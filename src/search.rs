use crate::core::{evoseed_keys, evoseed_seeds, EvoParams, SeedTrack};
use crate::index::SeedIndex;
use crate::scoring::{reject_score, weighted_score};
use anyhow::{bail, Result};
use std::collections::{BTreeMap, BTreeSet};

const MAX_POSITION_PAIRS_PER_KEY: usize = 64;

#[derive(Clone, Debug, PartialEq)]
pub struct CandidateScore {
    pub id: String,
    pub score: f64,
    pub shared: usize,
    pub chain_score: f64,
    pub weighted_jaccard: f64,
    pub raw_chain_score: f64,
    pub class_chain_score: f64,
}

#[derive(Clone, Debug)]
struct SeedHit {
    key: u128,
    track: SeedTrack,
    q1: usize,
    q2: usize,
    t1: usize,
    t2: usize,
    weight: f64,
}

impl SeedHit {
    fn q_center2(&self) -> i64 {
        (self.q1 + self.q2) as i64
    }

    fn t_center2(&self) -> i64 {
        (self.t1 + self.t2) as i64
    }

    fn diagonal2(&self) -> i64 {
        self.q_center2() - self.t_center2()
    }
}

#[derive(Default)]
struct Candidate {
    shared: BTreeSet<u128>,
    hits: Vec<SeedHit>,
    hit_counts: BTreeMap<u128, usize>,
}

struct Fenwick {
    tree: Vec<f64>,
}

impl Fenwick {
    fn new(n: usize) -> Self {
        Self {
            tree: vec![0.0; n + 1],
        }
    }

    /// 返回前 `end` 个元素的最大值，`end=0` 表示空前缀。
    fn prefix_max(&self, end: usize) -> f64 {
        let mut i = end;
        let mut best = 0.0_f64;
        while i > 0 {
            best = best.max(self.tree[i]);
            i -= i & (!i + 1);
        }
        best
    }

    fn update(&mut self, index: usize, value: f64) {
        let mut i = index;
        while i < self.tree.len() {
            self.tree[i] = self.tree[i].max(value);
            i += i & (!i + 1);
        }
    }
}

fn chain_on_offset(hits: &[SeedHit], track: SeedTrack, band: usize, offset: i64) -> f64 {
    let width = 2 * band as i64;
    let mut buckets: BTreeMap<i64, Vec<SeedHit>> = BTreeMap::new();
    for hit in hits.iter().filter(|hit| hit.track == track) {
        let bucket = (hit.diagonal2() + offset).div_euclid(width);
        buckets.entry(bucket).or_default().push(hit.clone());
    }
    let mut global_best = 0.0_f64;
    for mut bucket in buckets.into_values() {
        bucket.sort_by_key(|hit| (hit.q_center2(), hit.t_center2(), hit.key));
        let mut target_centers: Vec<i64> = bucket.iter().map(SeedHit::t_center2).collect();
        target_centers.sort_unstable();
        target_centers.dedup();
        let mut tree = Fenwick::new(target_centers.len());
        let mut start = 0;
        while start < bucket.len() {
            let q = bucket[start].q_center2();
            let mut end = start + 1;
            while end < bucket.len() && bucket[end].q_center2() == q {
                end += 1;
            }
            let mut updates = Vec::new();
            for hit in &bucket[start..end] {
                let t_index = target_centers.binary_search(&hit.t_center2()).unwrap();
                // 只查询严格更小的目标中心；同一查询中心稍后统一更新。
                let value = tree.prefix_max(t_index) + hit.weight;
                updates.push((t_index + 1, value));
                global_best = global_best.max(value);
            }
            for (index, value) in updates {
                tree.update(index, value);
            }
            start = end;
        }
    }
    global_best
}

/// 用两套错开的对角线分桶减轻硬桶边界影响，并在桶内求单调递增最优链。
fn best_chain(hits: &[SeedHit], track: SeedTrack, band: usize) -> f64 {
    chain_on_offset(hits, track, band, 0).max(chain_on_offset(hits, track, band, band as i64))
}

fn validate_min_shared(min_shared: usize) -> Result<()> {
    if min_shared == 0 {
        bail!("min_shared 必须 >= 1；零共享候选不属于倒排索引召回范围");
    }
    Ok(())
}

fn check_rejection<'a>(
    negative: &'a BTreeSet<u128>,
    db_keys: Option<&'a BTreeMap<String, BTreeSet<u128>>>,
) -> Result<Option<&'a BTreeMap<String, BTreeSet<u128>>>> {
    if !negative.is_empty() && db_keys.is_none() {
        bail!("负向排斥需要完整候选键集 db_keys");
    }
    Ok(db_keys)
}

/// 双通道、IDF 归一化的位置共线链检索。
///
/// 最终分数为 `最优链权重 / 查询与目标的加权键集并集`。每个唯一键在单个
/// 候选中的全部位置命中分摊一份 IDF 权重，限制重复序列刷分。
pub fn search(
    query: &str,
    index: &SeedIndex,
    params: &EvoParams,
    min_shared: usize,
    negative_keys: &BTreeSet<u128>,
    max_reject: usize,
    db_keys: Option<&BTreeMap<String, BTreeSet<u128>>>,
) -> Result<Vec<CandidateScore>> {
    validate_min_shared(min_shared)?;
    if params != &index.meta.params {
        bail!("查询参数与索引元数据不一致");
    }
    let db_keys = check_rejection(negative_keys, db_keys)?;
    let query_seeds = evoseed_seeds(query, params)?;
    let query_keys: BTreeSet<u128> = query_seeds.iter().map(|seed| seed.key).collect();
    let query_weight: f64 = query_keys.iter().map(|key| index.idf(key)).sum();
    let mut candidates: BTreeMap<String, Candidate> = BTreeMap::new();
    for seed in query_seeds {
        let Some(postings) = index.postings.get(&seed.key) else {
            continue;
        };
        for posting in postings {
            let candidate = candidates.entry(posting.id.clone()).or_default();
            candidate.shared.insert(seed.key);
            let count = candidate.hit_counts.get(&seed.key).copied().unwrap_or(0);
            if count < MAX_POSITION_PAIRS_PER_KEY {
                candidate.hits.push(SeedHit {
                    key: seed.key,
                    track: seed.track,
                    q1: seed.p1,
                    q2: seed.p2,
                    t1: posting.p1,
                    t2: posting.p2,
                    weight: 0.0,
                });
                candidate.hit_counts.insert(seed.key, count + 1);
            }
        }
    }

    let mut hits = Vec::new();
    for (id, candidate) in candidates {
        if candidate.shared.len() < min_shared {
            continue;
        }
        if !negative_keys.is_empty()
            && db_keys
                .and_then(|all| all.get(&id))
                .map(|keys| reject_score(keys, negative_keys) > max_reject)
                .unwrap_or(true)
        {
            continue;
        }
        let overlap: f64 = candidate.shared.iter().map(|key| index.idf(key)).sum();
        let target_weight = index
            .sequence_stats
            .get(&id)
            .map(|stat| stat.idf_weight)
            .unwrap_or(0.0);
        let union = query_weight + target_weight - overlap;
        let weighted_jaccard = if union > 0.0 { overlap / union } else { 0.0 };
        let Candidate {
            shared,
            hits: mut chain_hits,
            hit_counts,
        } = candidate;
        for hit in &mut chain_hits {
            hit.weight = index.idf(&hit.key) / hit_counts[&hit.key] as f64;
        }
        let raw_chain_score = best_chain(&chain_hits, SeedTrack::Raw, params.chain_band);
        let class_chain_score = best_chain(&chain_hits, SeedTrack::Class, params.chain_band);
        let chain_score = (raw_chain_score + class_chain_score).min(overlap);
        let score = if union > 0.0 {
            chain_score / union
        } else {
            0.0
        };
        hits.push(CandidateScore {
            id,
            score,
            shared: shared.len(),
            chain_score,
            weighted_jaccard,
            raw_chain_score,
            class_chain_score,
        });
    }
    hits.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then_with(|| b.raw_chain_score.total_cmp(&a.raw_chain_score))
            .then_with(|| b.weighted_jaccard.total_cmp(&a.weighted_jaccard))
            .then_with(|| b.shared.cmp(&a.shared))
            .then_with(|| a.id.cmp(&b.id))
    });
    Ok(hits)
}

/// 旧式唯一共享键计数，保留作消融对照；默认 CLI 不再使用。
pub fn search_count(
    query: &str,
    index: &SeedIndex,
    params: &EvoParams,
    min_shared: usize,
) -> Result<Vec<(String, usize)>> {
    validate_min_shared(min_shared)?;
    if params != &index.meta.params {
        bail!("查询参数与索引元数据不一致");
    }
    let query_keys = evoseed_keys(query, params)?;
    let mut votes: BTreeMap<String, usize> = BTreeMap::new();
    for key in query_keys {
        if let Some(postings) = index.postings.get(&key) {
            for id in postings.iter().map(|hit| &hit.id).collect::<BTreeSet<_>>() {
                *votes.entry(id.clone()).or_default() += 1;
            }
        }
    }
    let mut hits: Vec<_> = votes
        .into_iter()
        .filter(|(_, count)| *count >= min_shared)
        .collect();
    hits.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    Ok(hits)
}

pub fn search_weighted_index(
    query: &str,
    index: &SeedIndex,
    params: &EvoParams,
    weights: &BTreeMap<u128, f64>,
    min_score: Option<f64>,
    min_shared: usize,
) -> Result<Vec<(String, f64)>> {
    validate_min_shared(min_shared)?;
    if params != &index.meta.params {
        bail!("查询参数与索引元数据不一致");
    }
    if weights.values().any(|w| !w.is_finite())
        || min_score.map(|x| !x.is_finite()).unwrap_or(false)
    {
        bail!("权重和 min_score 必须为有限数值");
    }
    let query_keys = evoseed_keys(query, params)?;
    let mut votes: BTreeMap<String, (usize, f64)> = BTreeMap::new();
    for key in query_keys {
        let weight = weights.get(&key).copied().unwrap_or(0.0);
        if let Some(postings) = index.postings.get(&key) {
            for id in postings.iter().map(|hit| &hit.id).collect::<BTreeSet<_>>() {
                let entry = votes.entry(id.clone()).or_insert((0, 0.0));
                entry.0 += 1;
                entry.1 += weight;
            }
        }
    }
    let mut hits: Vec<_> = votes
        .into_iter()
        .filter(|(_, (shared, score))| {
            *shared >= min_shared
                && min_score
                    .map(|threshold| *score >= threshold)
                    .unwrap_or(true)
        })
        .map(|(id, (_, score))| (id, score))
        .collect();
    hits.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    Ok(hits)
}

pub fn search_with_keys(
    query: &str,
    db_keys: &BTreeMap<String, BTreeSet<u128>>,
    params: &EvoParams,
    min_shared: usize,
    negative_keys: &BTreeSet<u128>,
    max_reject: usize,
) -> Result<Vec<(String, usize)>> {
    validate_min_shared(min_shared)?;
    let query_keys = evoseed_keys(query, params)?;
    let mut hits: Vec<_> = db_keys
        .iter()
        .filter_map(|(id, keys)| {
            let shared = query_keys.intersection(keys).count();
            if shared < min_shared
                || (!negative_keys.is_empty() && reject_score(keys, negative_keys) > max_reject)
            {
                return None;
            }
            Some((id.clone(), shared))
        })
        .collect();
    hits.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    Ok(hits)
}

pub fn search_weighted_with_keys(
    query: &str,
    db_keys: &BTreeMap<String, BTreeSet<u128>>,
    params: &EvoParams,
    weights: &BTreeMap<u128, f64>,
    min_score: Option<f64>,
    min_shared: usize,
) -> Result<Vec<(String, f64)>> {
    validate_min_shared(min_shared)?;
    if weights.values().any(|w| !w.is_finite())
        || min_score.map(|x| !x.is_finite()).unwrap_or(false)
    {
        bail!("权重和 min_score 必须为有限数值");
    }
    let query_keys = evoseed_keys(query, params)?;
    let mut hits = Vec::new();
    for (id, keys) in db_keys {
        let shared = query_keys.intersection(keys).count();
        if shared < min_shared {
            continue;
        }
        let score = weighted_score(&query_keys, keys, weights, 0.0)?;
        if min_score
            .map(|threshold| score < threshold)
            .unwrap_or(false)
        {
            continue;
        }
        hits.push((id.clone(), score));
    }
    hits.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    Ok(hits)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::database_keys;

    const LONG: &str = "MVKAVLAVLAKGGTIAAQLSPEHMPVRMRMLDNNEAAELKALGVPFGRIRLPGGQVAEPIKVDANGVITDAQSSSLIFRGFERILQHQDPRDAAFLTQRICGVCPLSHGLTATNALDELYGVAEHVPKDGYKVTYAPQISDEEVKQMLEDLKRAG";
    const OTHER: &str = "MTNKKAAAWLVGAIVGGLGMAEAKLDDKAEAKKAAKEAKKAAELKALGVPFGRIRLPGGQVAEPIKVDANGVITDAQSSSLIFRGFERILQHQDPRDAAWQISNEEVKQMLEDLKRAG";

    fn setup() -> (BTreeMap<String, String>, SeedIndex, EvoParams) {
        let seq = "ACDEFGHIKLMNPQRSTVWY".repeat(8);
        let db = BTreeMap::from([("a".into(), seq.clone()), ("b".into(), seq)]);
        let params = EvoParams::default();
        let index = SeedIndex::build(&db, params.clone()).unwrap();
        (db, index, params)
    }

    fn synthetic_hit(key: u128, q: usize, t: usize) -> SeedHit {
        SeedHit {
            key,
            track: SeedTrack::Raw,
            q1: q,
            q2: q + 4,
            t1: t,
            t2: t + 4,
            weight: 1.0,
        }
    }

    #[test]
    fn collinear_chain_beats_scattered_hits() {
        let collinear = vec![
            synthetic_hit(1, 0, 0),
            synthetic_hit(2, 20, 20),
            synthetic_hit(3, 40, 40),
        ];
        let scattered = vec![
            synthetic_hit(1, 0, 40),
            synthetic_hit(2, 20, 20),
            synthetic_hit(3, 40, 0),
        ];
        assert_eq!(best_chain(&collinear, SeedTrack::Raw, 10), 3.0);
        assert_eq!(best_chain(&scattered, SeedTrack::Raw, 10), 1.0);
    }

    #[test]
    fn dual_chain_search_prefers_self_and_exposes_track_scores() {
        let db = BTreeMap::from([("self".into(), LONG.into()), ("other".into(), OTHER.into())]);
        let params = EvoParams::default();
        let index = SeedIndex::build(&db, params.clone()).unwrap();
        let hits = search(LONG, &index, &params, 1, &BTreeSet::new(), 0, None).unwrap();
        assert_eq!(hits[0].id, "self");
        assert!(hits[0].score > hits[1].score);
        assert!(hits[0].raw_chain_score > 0.0);
        assert!(hits[0].class_chain_score > 0.0);
    }

    #[test]
    fn zero_weight_candidates_and_zero_min_shared_are_handled() {
        let (db, index, params) = setup();
        let keys = database_keys(&db, &params).unwrap();
        let query = db["a"].clone();
        let weights = BTreeMap::new();
        let indexed =
            search_weighted_index(&query, &index, &params, &weights, Some(0.0), 1).unwrap();
        let brute =
            search_weighted_with_keys(&query, &keys, &params, &weights, Some(0.0), 1).unwrap();
        assert_eq!(indexed, brute);
        assert_eq!(indexed.len(), 2);
        assert!(search(&query, &index, &params, 0, &BTreeSet::new(), 0, None).is_err());
        assert!(search_with_keys(&query, &keys, &params, 0, &BTreeSet::new(), 0).is_err());
        let too_many = evoseed_keys(&query, &params).unwrap().len() + 1;
        let indexed =
            search_weighted_index(&query, &index, &params, &weights, Some(0.0), too_many).unwrap();
        let brute =
            search_weighted_with_keys(&query, &keys, &params, &weights, Some(0.0), too_many)
                .unwrap();
        assert!(indexed.is_empty());
        assert_eq!(indexed, brute);
        assert!(
            search_weighted_index(&query, &index, &params, &weights, Some(f64::NAN), 1).is_err()
        );
        assert!(
            search_weighted_with_keys(&query, &keys, &params, &weights, Some(f64::NAN), 1).is_err()
        );
    }

    #[test]
    fn count_index_and_brute_force_agree() {
        let (db, index, params) = setup();
        let keys = database_keys(&db, &params).unwrap();
        let query = db["a"].clone();
        let a = search_count(&query, &index, &params, 1).unwrap();
        let b = search_with_keys(&query, &keys, &params, 1, &BTreeSet::new(), 0).unwrap();
        assert_eq!(a, b);
        let mut mismatched = params.clone();
        mismatched.kc += 1;
        assert!(search(&query, &index, &mismatched, 1, &BTreeSet::new(), 0, None).is_err());
        assert!(
            search_weighted_index(&query, &index, &mismatched, &BTreeMap::new(), None, 1).is_err()
        );
    }
}
