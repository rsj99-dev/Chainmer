use crate::core::{evoseed_keys, EvoParams};
use crate::index::SeedIndex;
use crate::scoring::{reject_score, weighted_score};
use anyhow::{bail, Result};
use std::collections::{BTreeMap, BTreeSet};

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

pub fn search(
    query: &str,
    index: &SeedIndex,
    params: &EvoParams,
    min_shared: usize,
    negative_keys: &BTreeSet<u128>,
    max_reject: usize,
    db_keys: Option<&BTreeMap<String, BTreeSet<u128>>>,
) -> Result<Vec<(String, usize)>> {
    validate_min_shared(min_shared)?;
    if params != &index.meta.params {
        bail!("查询参数与索引元数据不一致");
    }
    let db_keys = check_rejection(negative_keys, db_keys)?;
    let query_keys = evoseed_keys(query, params)?;
    let mut votes: BTreeMap<String, usize> = BTreeMap::new();
    for key in query_keys {
        if let Some(ids) = index.postings.get(&key) {
            for id in ids {
                *votes.entry(id.clone()).or_default() += 1;
            }
        }
    }
    let mut hits: Vec<_> = votes
        .into_iter()
        .filter(|(id, count)| {
            if *count < min_shared {
                return false;
            }
            if negative_keys.is_empty() {
                return true;
            }
            db_keys
                .and_then(|all| all.get(id))
                .map(|keys| reject_score(keys, negative_keys) <= max_reject)
                .unwrap_or(false)
        })
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
        if let Some(ids) = index.postings.get(&key) {
            for id in ids {
                // 即使权重为零也建立候选，保持与穷举加权检索一致。
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

    fn setup() -> (BTreeMap<String, String>, SeedIndex, EvoParams) {
        let seq = "ACDEFGHIKLMNPQRSTVWY".repeat(8);
        let db = BTreeMap::from([("a".into(), seq.clone()), ("b".into(), seq)]);
        let params = EvoParams::default();
        let index = SeedIndex::build(&db, params.clone()).unwrap();
        (db, index, params)
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
    fn index_and_brute_force_agree() {
        let (db, index, params) = setup();
        let keys = database_keys(&db, &params).unwrap();
        let query = db["a"].clone();
        let a = search(&query, &index, &params, 1, &BTreeSet::new(), 0, None).unwrap();
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
