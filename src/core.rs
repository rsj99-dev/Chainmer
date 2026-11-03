use crate::seeds::{clean_protein, encode, hash_kmer, to_class_string, SecPolicy};
use anyhow::{bail, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
pub struct EvoParams {
    pub kc: usize,
    pub gmin: usize,
    pub gmax: usize,
    pub delta: usize,
    pub use_class: bool,
    pub gap_filter: bool,
    pub gap_bucket: bool,
    pub w1: usize,
    pub max_pairs_per_anchor: usize,
    pub sec_policy: String,
}

impl Default for EvoParams {
    fn default() -> Self {
        Self {
            kc: 5,
            gmin: 30,
            gmax: 300,
            delta: 15,
            use_class: true,
            gap_filter: true,
            gap_bucket: true,
            w1: 6,
            max_pairs_per_anchor: 12,
            sec_policy: "C".into(),
        }
    }
}

impl EvoParams {
    pub fn validate(&self) -> Result<()> {
        if self.kc == 0 || self.w1 == 0 || self.delta == 0 || self.max_pairs_per_anchor == 0 {
            bail!("kc、w1、delta 和 max_pairs_per_anchor 必须为正数");
        }
        if self.gmin > self.gmax {
            bail!("gmin 不能大于 gmax");
        }
        if !matches!(self.sec_policy.as_str(), "C" | "X") {
            bail!("sec_policy 仅支持 C 或 X");
        }
        Ok(())
    }

    pub fn set_use_gap(&mut self, enabled: bool) {
        self.gap_filter = enabled;
        self.gap_bucket = enabled;
    }

    fn sec_policy(&self) -> SecPolicy {
        if self.sec_policy == "C" {
            SecPolicy::C
        } else {
            SecPolicy::X
        }
    }
}

pub fn anchors(seq: &str, params: &EvoParams) -> Result<Vec<(usize, u128)>> {
    params.validate()?;
    let cleaned = clean_protein(seq, params.sec_policy());
    let transformed = if params.use_class {
        to_class_string(&cleaned)
    } else {
        cleaned
    };
    let codes = encode(&transformed);
    if codes.len() < params.kc {
        return Ok(Vec::new());
    }
    let km: Vec<(usize, u128)> = (0..=codes.len() - params.kc)
        .map(|i| {
            let part: Vec<u64> = codes[i..i + params.kc].iter().map(|x| *x as u64).collect();
            (i, hash_kmer(&part, 0xABCDEF))
        })
        .collect();
    if km.len() < params.w1 {
        return Ok(Vec::new());
    }

    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for i in 0..=km.len() - params.w1 {
        let best = *km[i..i + params.w1].iter().min_by_key(|x| x.1).unwrap();
        // 保留旧实现的确定性去重规则，确保 Rust/Python 可逐项核对。
        if seen.insert((best.1, best.0 / params.w1)) {
            out.push(best);
        }
    }
    Ok(out)
}

pub fn evoseed_keys(seq: &str, params: &EvoParams) -> Result<BTreeSet<u128>> {
    let anchors = anchors(seq, params)?;
    let mut keys = BTreeSet::new();
    for i in 0..anchors.len() {
        let (p1, h1) = anchors[i];
        let mut count = 0;
        for &(p2, h2) in anchors.iter().skip(i + 1) {
            if count >= params.max_pairs_per_anchor {
                break;
            }
            let d = p2 - p1;
            if params.gap_filter {
                if d < params.gmin {
                    continue;
                }
                if d > params.gmax {
                    break;
                }
            }
            let bucket = if params.gap_bucket {
                d / params.delta
            } else {
                0
            };
            // 与旧 API 一致：滚动累加时只保留各锚点哈希的低 64 bit。
            keys.insert(hash_kmer(&[h1 as u64, h2 as u64, bucket as u64], 0xE05EED));
            count += 1;
        }
    }
    Ok(keys)
}

#[cfg(test)]
mod tests {
    use super::*;

    const LONG: &str = "MVKAVLAVLAKGGTIAAQLSPEHMPVRMRMLDNNEAAELKALGVPFGRIRLPGGQVAEPIKVDANGVITDAQSSSLIFRGFERILQHQDPRDAAFLTQRICGVCPLSHGLTATNALDELYGVAEHVPKDGYKVTYAPQISDEEVKQMLEDLKRAG";

    #[test]
    fn gap_modes_and_pair_limit_are_real() {
        let mut p = EvoParams::default();
        let full = evoseed_keys(LONG, &p).unwrap();
        p.gap_bucket = false;
        let filter_only = evoseed_keys(LONG, &p).unwrap();
        p.set_use_gap(false);
        let no_gap = evoseed_keys(LONG, &p).unwrap();
        assert_ne!(full, filter_only);
        assert_ne!(filter_only, no_gap);
        p.max_pairs_per_anchor = 0;
        assert!(evoseed_keys(LONG, &p).is_err());
        p.max_pairs_per_anchor = 1;
        assert!(evoseed_keys(LONG, &p).unwrap().len() < no_gap.len());
    }

    #[test]
    fn invalid_parameters_fail_early() {
        let mut p = EvoParams {
            delta: 0,
            ..EvoParams::default()
        };
        assert!(evoseed_keys(LONG, &p).is_err());
        p = EvoParams::default();
        p.w1 = 0;
        assert!(evoseed_keys(LONG, &p).is_err());
        p = EvoParams::default();
        p.gmin = 20;
        p.gmax = 10;
        assert!(evoseed_keys(LONG, &p).is_err());
    }
}
