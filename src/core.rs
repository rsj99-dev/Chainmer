use crate::seeds::{clean_protein, encode, hash_kmer, spaced_kmers_stride, to_class_string, SecPolicy};
use anyhow::{bail, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SeedTrack {
    Raw,
    Class,
}

impl SeedTrack {
    fn salt(self) -> u64 {
        match self {
            Self::Raw => 0xE05EED,
            Self::Class => 0xE05EED ^ 0x000C_1A55_7AAC,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
pub struct PositionedSeed {
    pub key: u128,
    pub p1: usize,
    pub p2: usize,
    pub track: SeedTrack,
}

/// 密度混合通道：两个 10-11 位间距掩码（同 kmer_diamond 口径）。
// 常量已上移至 ChainmerParams 之前。

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
pub struct ChainmerParams {
    pub kc: usize,
    pub class_kc: usize,
    pub gmin: usize,
    pub gmax: usize,
    pub delta: usize,
    pub use_class: bool,
    pub dual_track: bool,
    pub gap_filter: bool,
    pub gap_bucket: bool,
    pub w1: usize,
    pub max_pairs_per_anchor: usize,
    pub chain_band: usize,
    pub sec_policy: String,
    /// 密度混合通道掩码；空 = 关闭（v1.1.0 起默认开启）。
    pub dense_masks: Vec<String>,
    /// 间距掩码步长（1 = 全位置采样；预分析实测 stride=1 达标）。
    pub dense_stride: usize,
    /// 是否叠加替换等价类轨道的间距键。
    pub dense_class: bool,
}

impl Default for ChainmerParams {
    fn default() -> Self {
        Self {
            kc: 5,
            class_kc: 8,
            gmin: 30,
            gmax: 300,
            delta: 5,
            use_class: true,
            dual_track: true,
            gap_filter: true,
            gap_bucket: true,
            w1: 6,
            max_pairs_per_anchor: 12,
            chain_band: 30,
            sec_policy: "C".into(),
            dense_masks: DENSE_MASKS.iter().map(|s| s.to_string()).collect(),
            dense_stride: 1,
            dense_class: true,
        }
    }
}

impl ChainmerParams {
    pub fn validate(&self) -> Result<()> {
        if self.kc == 0
            || self.class_kc == 0
            || self.w1 == 0
            || self.delta == 0
            || self.max_pairs_per_anchor == 0
            || self.chain_band == 0
        {
            bail!("kc、class_kc、w1、delta、max_pairs_per_anchor 和 chain_band 必须为正数");
        }
        if self.gmin > self.gmax {
            bail!("gmin 不能大于 gmax");
        }
        if !matches!(self.sec_policy.as_str(), "C" | "X") {
            bail!("sec_policy 仅支持 C 或 X");
        }
        if !self.dense_masks.is_empty() {
            if self.dense_stride == 0 {
                bail!("dense_stride 必须为正数");
            }
            for mask in &self.dense_masks {
                if mask.is_empty()
                    || !mask.bytes().all(|b| b == b'0' || b == b'1')
                    || !mask.contains('1')
                {
                    bail!("dense_masks 必须是含 1 的非空 0/1 字符串");
                }
            }
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

fn track_anchors_impl(
    seq: &str,
    kc: usize,
    w1: usize,
    track: SeedTrack,
    sec_policy: SecPolicy,
    mask_unknown: bool,
) -> Vec<(usize, u128)> {
    let cleaned = clean_protein(seq, sec_policy);
    let valid: Vec<bool> = cleaned.bytes().map(|residue| residue != b'X').collect();
    let transformed = match track {
        SeedTrack::Raw => cleaned,
        SeedTrack::Class => to_class_string(&cleaned),
    };
    let codes = encode(&transformed);
    if codes.len() < kc || w1 == 0 {
        return Vec::new();
    }
    let km: Vec<(usize, u128)> = (0..=codes.len() - kc)
        .filter(|i| !mask_unknown || valid[*i..*i + kc].iter().all(|value| *value))
        .map(|i| {
            let part: Vec<u64> = codes[i..i + kc].iter().map(|x| *x as u64).collect();
            (i, hash_kmer(&part, 0xABCDEF))
        })
        .collect();
    if km.len() < w1 {
        return Vec::new();
    }

    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for i in 0..=km.len() - w1 {
        let best = *km[i..i + w1].iter().min_by_key(|x| x.1).unwrap();
        if seen.insert((best.1, best.0 / w1)) {
            out.push(best);
        }
    }
    out
}

fn track_anchors(
    seq: &str,
    kc: usize,
    w1: usize,
    track: SeedTrack,
    sec_policy: SecPolicy,
) -> Vec<(usize, u128)> {
    track_anchors_impl(seq, kc, w1, track, sec_policy, false)
}

/// 兼容旧的单通道锚点接口；新检索使用 `chainmer_seeds` 的双通道位置种子。
pub fn anchors(seq: &str, params: &ChainmerParams) -> Result<Vec<(usize, u128)>> {
    params.validate()?;
    let track = if params.use_class {
        SeedTrack::Class
    } else {
        SeedTrack::Raw
    };
    Ok(track_anchors(
        seq,
        params.kc,
        params.w1,
        track,
        params.sec_policy(),
    ))
}

/// 新主通道默认屏蔽所有覆盖 X/非标准残基的锚点窗口。
pub fn chainmer_track_seeds(
    seq: &str,
    params: &ChainmerParams,
    track: SeedTrack,
) -> Result<BTreeSet<PositionedSeed>> {
    chainmer_track_seeds_impl(seq, params, track, true)
}

fn chainmer_track_seeds_impl(
    seq: &str,
    params: &ChainmerParams,
    track: SeedTrack,
    mask_unknown: bool,
) -> Result<BTreeSet<PositionedSeed>> {
    params.validate()?;
    let kc = match track {
        SeedTrack::Raw => params.kc,
        SeedTrack::Class => params.class_kc,
    };
    let anchors = track_anchors_impl(seq, kc, params.w1, track, params.sec_policy(), mask_unknown);
    let mut seeds = BTreeSet::new();
    for (i, &(p1, h1)) in anchors.iter().enumerate() {
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
            seeds.insert(PositionedSeed {
                key: hash_kmer(&[h1 as u64, h2 as u64, bucket as u64], track.salt()),
                p1,
                p2,
                track,
            });
            count += 1;
        }
    }
    Ok(seeds)
}

/// 生成 raw/class 双通道位置种子；`dual_track=false` 时按 `use_class` 选择单通道。
pub fn chainmer_seeds(seq: &str, params: &ChainmerParams) -> Result<BTreeSet<PositionedSeed>> {
    params.validate()?;
    let mut out = BTreeSet::new();
    if params.dual_track {
        out.extend(chainmer_track_seeds(seq, params, SeedTrack::Raw)?);
        out.extend(chainmer_track_seeds(seq, params, SeedTrack::Class)?);
    } else {
        let track = if params.use_class {
            SeedTrack::Class
        } else {
            SeedTrack::Raw
        };
        out.extend(chainmer_track_seeds(seq, params, track)?);
    }
    Ok(out)
}

/// 密度混合通道：两个 10-11 位间距掩码（同 kmer_diamond 口径），默认步长 2。
pub const DENSE_MASKS: [&str; 2] = ["111101110111", "111011010010111"];

/// 间距 k-mer 加密键（raw 残基或替换等价类轨道）。
pub fn dense_spaced_keys_track(
    seq: &str,
    masks: &[&str],
    stride: usize,
    class_track: bool,
) -> Result<BTreeSet<u128>> {
    let cleaned = clean_protein(seq, SecPolicy::C);
    let transformed = if class_track {
        to_class_string(&cleaned)
    } else {
        cleaned
    };
    let salt = if class_track { 0x0C1A55 } else { 0x0DE25E };
    let mut out = BTreeSet::new();
    for mask in masks {
        let k = mask.bytes().filter(|b| *b == b'1').count();
        for (_, h) in spaced_kmers_stride(&transformed, mask, k, stride)? {
            out.insert(hash_kmer(&[(h & u64::MAX as u128) as u64], salt));
        }
    }
    Ok(out)
}

pub fn dense_spaced_keys(
    seq: &str,
    masks: &[&str],
    stride: usize,
) -> Result<BTreeSet<u128>> {
    dense_spaced_keys_track(seq, masks, stride, false)
}

/// 仅间距加密键（不含锚点对键）；索引 v7 的 per-target dense 表用它。
pub fn chainmer_dense_keys(seq: &str, params: &ChainmerParams) -> Result<BTreeSet<u128>> {
    if params.dense_masks.is_empty() {
        return Ok(BTreeSet::new());
    }
    let masks: Vec<&str> = params.dense_masks.iter().map(|s| s.as_str()).collect();
    let stride = params.dense_stride.max(1);
    let mut out = dense_spaced_keys_track(seq, &masks, stride, false)?;
    if params.dense_class {
        out.extend(dense_spaced_keys_track(seq, &masks, stride, true)?);
    }
    Ok(out)
}

/// 双通道锚点对键 ∪ 间距加密键（密度混合默认路径）。
pub fn chainmer_keys_dense(seq: &str, params: &ChainmerParams) -> Result<BTreeSet<u128>> {
    chainmer_keys_dense_full(seq, params, 1, true)
}

pub fn chainmer_keys(seq: &str, params: &ChainmerParams) -> Result<BTreeSet<u128>> {
    chainmer_keys_dense_full(
        seq,
        params,
        params.dense_stride.max(1),
        params.dense_class,
    )
}

/// 升序去重 Vec 版（免 BTreeSet 节点分配；检索/预计算热路径用）。
pub fn chainmer_keys_vec(seq: &str, params: &ChainmerParams) -> Result<Vec<u128>> {
    let mut v: Vec<u128> = chainmer_seeds(seq, params)?
        .into_iter()
        .map(|s| s.key)
        .collect();
    if !params.dense_masks.is_empty() {
        let masks: Vec<&str> = params.dense_masks.iter().map(|s| s.as_str()).collect();
        let stride = params.dense_stride.max(1);
        v.extend(dense_spaced_keys_track(seq, &masks, stride, false)?);
        if params.dense_class {
            v.extend(dense_spaced_keys_track(seq, &masks, stride, true)?);
        }
    }
    v.sort_unstable();
    v.dedup();
    Ok(v)
}

/// 双通道锚点对键 ∪ 间距加密键；stride/class 可显式覆盖（消融与环境变量调参用）。
pub fn chainmer_keys_dense_full(
    seq: &str,
    params: &ChainmerParams,
    stride: usize,
    class_dense: bool,
) -> Result<BTreeSet<u128>> {
    let mut out: BTreeSet<u128> = chainmer_seeds(seq, params)?
        .into_iter()
        .map(|seed| seed.key)
        .collect();
    if !params.dense_masks.is_empty() {
        let masks: Vec<&str> = params.dense_masks.iter().map(|s| s.as_str()).collect();
        out.extend(dense_spaced_keys_track(seq, &masks, stride, false)?);
        if class_dense {
            out.extend(dense_spaced_keys_track(seq, &masks, stride, true)?);
        }
    }
    Ok(out)
}

/// 仅锚点对键（v1.0.x 口径），保留作消融对照。
pub fn chainmer_pair_keys(seq: &str, params: &ChainmerParams) -> Result<BTreeSet<u128>> {
    Ok(chainmer_seeds(seq, params)?
        .into_iter()
        .map(|seed| seed.key)
        .collect())
}

/// v0.4 单通道键生成，仅用于历史对照，不参与新索引和默认检索。
pub fn chainmer_keys_legacy(seq: &str, params: &ChainmerParams) -> Result<BTreeSet<u128>> {
    params.validate()?;
    let track = if params.use_class {
        SeedTrack::Class
    } else {
        SeedTrack::Raw
    };
    let anchors = track_anchors(seq, params.kc, params.w1, track, params.sec_policy());
    let mut keys = BTreeSet::new();
    for (i, &(p1, h1)) in anchors.iter().enumerate() {
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
    fn dual_track_has_disjoint_keys_and_real_positions() {
        let p = ChainmerParams::default();
        let raw = chainmer_track_seeds(LONG, &p, SeedTrack::Raw).unwrap();
        let class = chainmer_track_seeds(LONG, &p, SeedTrack::Class).unwrap();
        assert!(!raw.is_empty() && !class.is_empty());
        assert!(raw.iter().all(|s| s.track == SeedTrack::Raw && s.p2 > s.p1));
        assert!(class.iter().all(|s| s.track == SeedTrack::Class));
        assert!(raw
            .iter()
            .map(|s| s.key)
            .collect::<BTreeSet<_>>()
            .is_disjoint(&class.iter().map(|s| s.key).collect::<BTreeSet<_>>()));
        let both = chainmer_seeds(LONG, &p).unwrap();
        assert_eq!(both.len(), raw.len() + class.len());
    }

    #[test]
    fn gap_modes_pair_limit_and_class_length_are_real() {
        let mut p = ChainmerParams::default();
        let full = chainmer_keys(LONG, &p).unwrap();
        p.gap_bucket = false;
        let filter_only = chainmer_keys(LONG, &p).unwrap();
        p.set_use_gap(false);
        let no_gap = chainmer_keys(LONG, &p).unwrap();
        assert_ne!(full, filter_only);
        assert_ne!(filter_only, no_gap);
        p.max_pairs_per_anchor = 0;
        assert!(chainmer_keys(LONG, &p).is_err());
        p.max_pairs_per_anchor = 1;
        assert!(chainmer_keys(LONG, &p).unwrap().len() < no_gap.len());
        let mut p = ChainmerParams::default();
        let long_class = chainmer_track_seeds(LONG, &p, SeedTrack::Class).unwrap();
        p.class_kc = 5;
        assert_ne!(
            long_class,
            chainmer_track_seeds(LONG, &p, SeedTrack::Class).unwrap()
        );
    }

    #[test]
    fn unknown_residues_are_masked_only_in_the_current_channel() {
        let p = ChainmerParams::default();
        let standard_raw = chainmer_track_seeds(LONG, &p, SeedTrack::Raw).unwrap();
        let reference_raw = chainmer_track_seeds_impl(LONG, &p, SeedTrack::Raw, false).unwrap();
        assert_eq!(standard_raw, reference_raw);
        assert!(chainmer_seeds(&"X".repeat(400), &p).unwrap().is_empty());
        assert!(chainmer_seeds(&"BZJO*".repeat(80), &p).unwrap().is_empty());
        let with_x = format!("{}X{}", &LONG[..80], &LONG[81..]);
        let raw = track_anchors_impl(&with_x, p.kc, p.w1, SeedTrack::Raw, SecPolicy::C, true);
        let class = track_anchors_impl(
            &with_x,
            p.class_kc,
            p.w1,
            SeedTrack::Class,
            SecPolicy::C,
            true,
        );
        assert!(raw
            .iter()
            .all(|(position, _)| *position + p.kc <= 80 || *position >= 81));
        assert!(class
            .iter()
            .all(|(position, _)| *position + p.class_kc <= 80 || *position >= 81));
        let legacy = chainmer_keys_legacy(&with_x, &p).unwrap();
        assert!(!legacy.is_empty());
        let u_as_c = track_anchors_impl("ACDU", 2, 1, SeedTrack::Raw, SecPolicy::C, true);
        let u_as_x = track_anchors_impl("ACDU", 2, 1, SeedTrack::Raw, SecPolicy::X, true);
        assert_eq!(u_as_c.len(), 3);
        assert_eq!(u_as_x.len(), 2);
    }

    #[test]
    fn extreme_parameters_return_without_arithmetic_panics() {
        let mut p = ChainmerParams {
            kc: usize::MAX,
            class_kc: usize::MAX,
            ..ChainmerParams::default()
        };
        assert!(chainmer_seeds(LONG, &p).unwrap().is_empty());
        p = ChainmerParams {
            w1: usize::MAX,
            ..ChainmerParams::default()
        };
        assert!(chainmer_seeds(LONG, &p).unwrap().is_empty());
        p = ChainmerParams {
            gmin: usize::MAX,
            gmax: usize::MAX,
            delta: usize::MAX,
            max_pairs_per_anchor: usize::MAX,
            chain_band: usize::MAX,
            ..ChainmerParams::default()
        };
        p.validate().unwrap();
        assert!(chainmer_seeds(LONG, &p).unwrap().is_empty());
    }

    #[test]
    fn invalid_parameters_fail_early() {
        let mut p = ChainmerParams {
            delta: 0,
            ..ChainmerParams::default()
        };
        assert!(chainmer_keys(LONG, &p).is_err());
        p = ChainmerParams::default();
        p.w1 = 0;
        assert!(chainmer_keys(LONG, &p).is_err());
        p = ChainmerParams::default();
        p.gmin = 20;
        p.gmax = 10;
        assert!(chainmer_keys(LONG, &p).is_err());
        p = ChainmerParams::default();
        p.chain_band = 0;
        assert!(chainmer_keys(LONG, &p).is_err());
    }
}
