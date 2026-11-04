use crate::core::{chainmer_keys, chainmer_seeds, ChainmerParams, SeedTrack};
use crate::index::SeedIndex;
use crate::scoring::{reject_score, weighted_score};
use crate::seeds::hash_kmer;
use anyhow::{bail, Context, Result};
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
    pub position_hits_total: usize,
    pub position_hits_retained: usize,
    pub truncated_keys: usize,
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
    fn q_center2(&self) -> i128 {
        self.q1 as i128 + self.q2 as i128
    }

    fn t_center2(&self) -> i128 {
        self.t1 as i128 + self.t2 as i128
    }

    fn diagonal2(&self) -> i128 {
        self.q_center2() - self.t_center2()
    }
}

/// 四个端点必须在两条序列中沿相同方向严格前进；禁止同角色端点重复和交叉映射。
fn can_chain(previous: &SeedHit, current: &SeedHit) -> bool {
    previous.q1 < current.q1
        && previous.q2 < current.q2
        && previous.t1 < current.t1
        && previous.t2 < current.t2
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct HitOrder {
    rank: u128,
    q1: usize,
    q2: usize,
    t1: usize,
    t2: usize,
}

impl HitOrder {
    fn new(hit: &SeedHit) -> Self {
        let rank = hash_kmer(
            &[
                (hit.key >> 64) as u64,
                hit.key as u64,
                hit.q1 as u64,
                hit.q2 as u64,
                hit.t1 as u64,
                hit.t2 as u64,
            ],
            0xC1A1_9EED,
        );
        Self {
            rank,
            q1: hit.q1,
            q2: hit.q2,
            t1: hit.t1,
            t2: hit.t2,
        }
    }
}

/// 稳定哈希 bottom-k：保留集合与 posting 遍历顺序无关，且不偏向序列任一端。
#[derive(Default)]
struct RetainedHits {
    total: usize,
    ordered: BTreeMap<HitOrder, SeedHit>,
}

impl RetainedHits {
    fn offer(&mut self, hit: SeedHit) -> Result<()> {
        self.total = self
            .total
            .checked_add(1)
            .context("位置命中计数超出 usize 范围")?;
        let order = HitOrder::new(&hit);
        if self.ordered.len() < MAX_POSITION_PAIRS_PER_KEY {
            self.ordered.insert(order, hit);
        } else if self
            .ordered
            .last_key_value()
            .map(|(last, _)| order < *last)
            .unwrap_or(true)
        {
            self.ordered.pop_last();
            self.ordered.insert(order, hit);
        }
        Ok(())
    }
}

#[derive(Default)]
struct Candidate {
    shared: BTreeSet<u128>,
    hits: BTreeMap<u128, RetainedHits>,
}

fn chain_on_offset(hits: &[SeedHit], track: SeedTrack, band: usize, offset: i128) -> f64 {
    let width = 2 * band as i128;
    let mut buckets: BTreeMap<i128, Vec<SeedHit>> = BTreeMap::new();
    for hit in hits.iter().filter(|hit| hit.track == track) {
        let bucket = (hit.diagonal2() + offset).div_euclid(width);
        buckets.entry(bucket).or_default().push(hit.clone());
    }
    let mut global_best = 0.0_f64;
    for mut bucket in buckets.into_values() {
        bucket.sort_by_key(|hit| (hit.q1, hit.q2, hit.t1, hit.t2, hit.key));
        let mut dp = vec![0.0_f64; bucket.len()];
        for i in 0..bucket.len() {
            dp[i] = bucket[i].weight;
            for j in 0..i {
                if can_chain(&bucket[j], &bucket[i]) {
                    dp[i] = dp[i].max(dp[j] + bucket[i].weight);
                }
            }
            global_best = global_best.max(dp[i]);
        }
    }
    global_best
}

/// 两套对角线桶仅用于容差筛选；桶内以严格端点对应关系求最优链。
fn best_chain(hits: &[SeedHit], track: SeedTrack, band: usize) -> f64 {
    chain_on_offset(hits, track, band, 0).max(chain_on_offset(hits, track, band, band as i128))
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
/// 候选中的保留代表命中共同分摊一份 IDF 权重，限制重复序列刷分。
pub fn search(
    query: &str,
    index: &SeedIndex,
    params: &ChainmerParams,
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
    let query_seeds = chainmer_seeds(query, params)?;
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
            candidate.hits.entry(seed.key).or_default().offer(SeedHit {
                key: seed.key,
                track: seed.track,
                q1: seed.p1,
                q2: seed.p2,
                t1: posting.p1,
                t2: posting.p2,
                weight: 0.0,
            })?;
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
            hits: retained_by_key,
        } = candidate;
        let position_hits_total = retained_by_key.values().map(|hits| hits.total).sum();
        let truncated_keys = retained_by_key
            .values()
            .filter(|hits| hits.total > hits.ordered.len())
            .count();
        let mut chain_hits = Vec::new();
        for (key, retained) in retained_by_key {
            let retained_count = retained.ordered.len();
            let weight = index.idf(&key) / retained_count as f64;
            chain_hits.extend(retained.ordered.into_values().map(|mut hit| {
                hit.weight = weight;
                hit
            }));
        }
        let position_hits_retained = chain_hits.len();
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
            position_hits_total,
            position_hits_retained,
            truncated_keys,
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
    params: &ChainmerParams,
    min_shared: usize,
) -> Result<Vec<(String, usize)>> {
    validate_min_shared(min_shared)?;
    if params != &index.meta.params {
        bail!("查询参数与索引元数据不一致");
    }
    let query_keys = chainmer_keys(query, params)?;
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
    params: &ChainmerParams,
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
    let query_keys = chainmer_keys(query, params)?;
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
    params: &ChainmerParams,
    min_shared: usize,
    negative_keys: &BTreeSet<u128>,
    max_reject: usize,
) -> Result<Vec<(String, usize)>> {
    validate_min_shared(min_shared)?;
    let query_keys = chainmer_keys(query, params)?;
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
    params: &ChainmerParams,
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
    let query_keys = chainmer_keys(query, params)?;
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

    fn setup() -> (BTreeMap<String, String>, SeedIndex, ChainmerParams) {
        let seq = "ACDEFGHIKLMNPQRSTVWY".repeat(8);
        let db = BTreeMap::from([("a".into(), seq.clone()), ("b".into(), seq)]);
        let params = ChainmerParams::default();
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

    fn endpoint_hit(q1: usize, q2: usize, t1: usize, t2: usize) -> SeedHit {
        SeedHit {
            key: 1,
            track: SeedTrack::Raw,
            q1,
            q2,
            t1,
            t2,
            weight: 1.0,
        }
    }

    fn reference_chain_on_offset(
        hits: &[SeedHit],
        track: SeedTrack,
        band: usize,
        offset: i128,
    ) -> f64 {
        let width = 2 * band as i128;
        let mut buckets: BTreeMap<i128, Vec<SeedHit>> = BTreeMap::new();
        for hit in hits.iter().filter(|hit| hit.track == track) {
            let bucket = (hit.diagonal2() + offset).div_euclid(width);
            buckets.entry(bucket).or_default().push(hit.clone());
        }
        let mut global_best = 0.0_f64;
        for bucket in buckets.into_values() {
            assert!(bucket.len() <= 8);
            for mask in 1_usize..(1 << bucket.len()) {
                let mut selected: Vec<SeedHit> = (0..bucket.len())
                    .filter(|i| mask & (1 << i) != 0)
                    .map(|i| bucket[i].clone())
                    .collect();
                selected.sort_by_key(|hit| (hit.q1, hit.q2, hit.t1, hit.t2, hit.key));
                if selected
                    .windows(2)
                    .all(|pair| can_chain(&pair[0], &pair[1]))
                {
                    global_best = global_best.max(selected.iter().map(|hit| hit.weight).sum());
                }
            }
        }
        global_best
    }

    #[test]
    fn endpoint_geometry_rejects_crossed_and_repeated_anchors() {
        let a = endpoint_hit(0, 100, 0, 100);
        let crossed = endpoint_hit(25, 90, 36, 101);
        let compatible = endpoint_hit(25, 110, 36, 120);
        let hits = vec![a.clone(), crossed];
        assert_eq!(chain_on_offset(&hits, SeedTrack::Raw, 30, 30), 1.0);
        let hits = vec![a.clone(), compatible];
        assert_eq!(chain_on_offset(&hits, SeedTrack::Raw, 30, 30), 2.0);
        let same_role_shared = endpoint_hit(0, 120, 0, 120);
        assert!(!can_chain(&a, &same_role_shared));
        let touching = endpoint_hit(100, 120, 100, 120);
        assert!(can_chain(&a, &touching));
        let same_center = endpoint_hit(25, 75, 25, 75);
        assert!(!can_chain(&a, &same_center));
    }

    #[test]
    fn diagonal_bucket_boundaries_and_extreme_band_are_defined() {
        let a = endpoint_hit(0, 100, 0, 100);
        let inside = endpoint_hit(200, 300, 171, 271);
        assert_eq!(
            chain_on_offset(&[a.clone(), inside], SeedTrack::Raw, 30, 0),
            2.0
        );
        let boundary = endpoint_hit(200, 300, 170, 270);
        assert_eq!(
            chain_on_offset(&[a.clone(), boundary], SeedTrack::Raw, 30, 0),
            1.0
        );
        assert_eq!(best_chain(&[a], SeedTrack::Raw, usize::MAX), 1.0);
    }

    #[test]
    fn chain_dp_matches_independent_subset_oracle() {
        let mut state = 17_u64;
        for _ in 0..40 {
            let mut hits = Vec::new();
            for _ in 0..7 {
                let mut next = || {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    (state % 45) as usize
                };
                let q1 = next();
                let q2 = q1 + 1 + next() % 15;
                let t1 = next();
                let t2 = t1 + 1 + next() % 15;
                hits.push(endpoint_hit(q1, q2, t1, t2));
            }
            for band in [3, 10, 30] {
                for offset in [0, band as i128] {
                    assert_eq!(
                        chain_on_offset(&hits, SeedTrack::Raw, band, offset),
                        reference_chain_on_offset(&hits, SeedTrack::Raw, band, offset)
                    );
                }
            }
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
    fn position_retention_is_order_independent_and_reports_truncation() {
        let query = "ACDEFGHIKLMNPQRSTVWY".repeat(40);
        let db = BTreeMap::from([("a".into(), query.clone())]);
        let params = ChainmerParams {
            dual_track: false,
            use_class: false,
            ..ChainmerParams::default()
        };
        let index = SeedIndex::build(&db, params.clone()).unwrap();
        let baseline = search(&query, &index, &params, 1, &BTreeSet::new(), 0, None).unwrap();
        assert!(baseline.iter().any(|hit| {
            hit.position_hits_total > hit.position_hits_retained && hit.truncated_keys > 0
        }));
        assert!(baseline
            .iter()
            .all(|hit| hit.position_hits_retained <= 64 * hit.shared));
        let mut reversed = index;
        for postings in reversed.postings.values_mut() {
            postings.reverse();
        }
        let shuffled = search(&query, &reversed, &params, 1, &BTreeSet::new(), 0, None).unwrap();
        assert_eq!(shuffled, baseline);
    }

    #[test]
    fn representative_retention_can_preserve_a_late_optimal_chain() {
        let early: Vec<_> = (0..80)
            .map(|i| endpoint_hit(i * 10, i * 10 + 5, 1_000_000 - i * 20, 1_000_005 - i * 20))
            .collect();
        for start in (10_000..20_000).step_by(100) {
            let late = vec![
                endpoint_hit(start, start + 5, start, start + 5),
                endpoint_hit(start + 20, start + 25, start + 20, start + 25),
            ];
            let mut all = early.clone();
            all.extend(late.clone());
            let mut retained = RetainedHits::default();
            for hit in &all {
                retained.offer(hit.clone()).unwrap();
            }
            let selected: Vec<_> = retained.ordered.into_values().collect();
            if late.iter().all(|late_hit| {
                selected.iter().any(|hit| {
                    (hit.q1, hit.q2, hit.t1, hit.t2)
                        == (late_hit.q1, late_hit.q2, late_hit.t1, late_hit.t2)
                })
            }) {
                assert_eq!(best_chain(&all, SeedTrack::Raw, 30), 2.0);
                assert_eq!(best_chain(&selected, SeedTrack::Raw, 30), 2.0);
                return;
            }
        }
        panic!("代表保留在所有测试位置上都丢失了后段最优链");
    }

    #[test]
    fn dual_chain_search_prefers_self_and_exposes_track_scores() {
        let db = BTreeMap::from([("self".into(), LONG.into()), ("other".into(), OTHER.into())]);
        let params = ChainmerParams::default();
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
        let too_many = chainmer_keys(&query, &params).unwrap().len() + 1;
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
