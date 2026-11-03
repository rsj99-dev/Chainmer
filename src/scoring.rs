use anyhow::{bail, Result};
use std::collections::{BTreeMap, BTreeSet};

pub fn consistency_score(a: &BTreeSet<u128>, b: &BTreeSet<u128>) -> usize {
    a.intersection(b).count()
}

pub fn reject_score(candidate: &BTreeSet<u128>, negative: &BTreeSet<u128>) -> usize {
    candidate.intersection(negative).count()
}

pub fn build_negative_set(pair_intersections: &[BTreeSet<u128>]) -> BTreeSet<u128> {
    let mut freq: BTreeMap<u128, usize> = BTreeMap::new();
    for keys in pair_intersections {
        for key in keys {
            *freq.entry(*key).or_default() += 1;
        }
    }
    freq.into_iter()
        .filter_map(|(k, n)| (n >= 2).then_some(k))
        .collect()
}

pub fn specificity_weights(
    pair_intersections: &[BTreeSet<u128>],
    n_neg: usize,
) -> Result<BTreeMap<u128, f64>> {
    if n_neg == 0 {
        bail!("n_neg 必须为正数");
    }
    let mut freq: BTreeMap<u128, usize> = BTreeMap::new();
    for keys in pair_intersections {
        for key in keys {
            *freq.entry(*key).or_default() += 1;
        }
    }
    Ok(freq
        .into_iter()
        .map(|(k, n)| (k, ((n_neg + 1) as f64 / (n + 1) as f64).ln()))
        .collect())
}

pub fn contrastive_weights(
    positive: &[BTreeSet<u128>],
    negative: &[BTreeSet<u128>],
) -> BTreeMap<u128, f64> {
    let mut pos: BTreeMap<u128, usize> = BTreeMap::new();
    let mut neg: BTreeMap<u128, usize> = BTreeMap::new();
    for keys in positive {
        for key in keys {
            *pos.entry(*key).or_default() += 1;
        }
    }
    for keys in negative {
        for key in keys {
            *neg.entry(*key).or_default() += 1;
        }
    }
    let all: BTreeSet<u128> = pos.keys().chain(neg.keys()).copied().collect();
    all.into_iter()
        .map(|k| {
            let score = ((pos.get(&k).copied().unwrap_or(0) + 1) as f64
                / (neg.get(&k).copied().unwrap_or(0) + 1) as f64)
                .ln();
            (k, score)
        })
        .collect()
}

pub fn weighted_score(
    query: &BTreeSet<u128>,
    candidate: &BTreeSet<u128>,
    weights: &BTreeMap<u128, f64>,
    default_weight: f64,
) -> Result<f64> {
    if weights.values().any(|w| !w.is_finite()) || !default_weight.is_finite() {
        bail!("权重必须为有限数值");
    }
    Ok(query
        .intersection(candidate)
        .map(|k| weights.get(k).copied().unwrap_or(default_weight))
        .sum())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contrastive_direction_and_neutral_default() {
        let pos = vec![BTreeSet::from([1, 2])];
        let neg = vec![BTreeSet::from([2, 3])];
        let w = contrastive_weights(&pos, &neg);
        assert!(w[&1] > 0.0 && w[&3] < 0.0);
        assert_eq!(w[&2], 0.0);
        assert_eq!(
            weighted_score(&BTreeSet::from([9]), &BTreeSet::from([9]), &w, 0.0).unwrap(),
            0.0
        );
    }
}
