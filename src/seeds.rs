use anyhow::{bail, Result};
use std::collections::{BTreeMap, BTreeSet};

pub const AA20: &str = "ACDEFGHIKLMNPQRSTVWY";
pub const AMBIGUOUS: u8 = 20;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SecPolicy {
    C,
    X,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SyncmerMode {
    Closed,
    Open { offset: usize },
}

pub fn clean_protein(seq: &str, sec_policy: SecPolicy) -> String {
    seq.chars()
        .filter_map(|raw| {
            let mut ch = raw.to_ascii_uppercase();
            if matches!(ch, ' ' | '\t' | '\r' | '\n' | '-' | '.') {
                return None;
            }
            if ch == 'U' {
                ch = match sec_policy {
                    SecPolicy::C => 'C',
                    SecPolicy::X => 'X',
                };
            }
            if ch == '*' {
                ch = 'X';
            }
            Some(if AA20.contains(ch) { ch } else { 'X' })
        })
        .collect()
}

pub fn encode(seq: &str) -> Vec<u8> {
    seq.chars()
        .map(|ch| match ch {
            'A'..='Z' => AA20.find(ch).map(|x| x as u8).unwrap_or(AMBIGUOUS),
            'a'..='f' => 21 + (ch as u8 - b'a'),
            _ => AMBIGUOUS,
        })
        .collect()
}

pub fn class_map() -> BTreeMap<String, String> {
    let groups = [
        ("a", "AVLIMC"),
        ("b", "STNPG"),
        ("c", "DEQ"),
        ("d", "KRH"),
        ("e", "FYW"),
        ("f", "XBZJOU*"),
    ];
    let mut map = BTreeMap::new();
    for (class, residues) in groups {
        for residue in residues.chars() {
            map.insert(residue.to_string(), class.to_string());
        }
    }
    map
}

pub fn to_class_string(seq: &str) -> String {
    let map = class_map();
    seq.chars()
        .map(|ch| map.get(&ch.to_string()).map(|x| x.as_str()).unwrap_or("f"))
        .collect()
}

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E3779B97F4A7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}

/// 与 Python 原型逐项兼容的哈希。
///
/// 输出以 u128 存储，但 hi/lo 均来自同一个 64-bit 状态，故有效状态最多 64 bit。
/// 不得将其描述为具有 128-bit 熵或 128-bit 碰撞保证。
pub fn hash_kmer(codes: &[u64], salt: u64) -> u128 {
    let mut h = salt;
    for code in codes {
        h = (h << 5).wrapping_add(h).wrapping_add(*code).wrapping_add(1);
    }
    let lo = splitmix64(h ^ 0x9E3779B97F4A7C15);
    let hi = splitmix64(h ^ 0xD1B54A32D192ED03);
    ((hi as u128) << 64) | lo as u128
}

pub fn hash_bytes(bytes: &[u8], salt: u64) -> u128 {
    let codes: Vec<u64> = bytes.iter().map(|x| *x as u64).collect();
    hash_kmer(&codes, salt)
}

pub fn kmers(seq: &str, k: usize) -> Result<Vec<(usize, u128)>> {
    if k == 0 {
        bail!("k 必须为正数");
    }
    let codes = encode(seq);
    if codes.len() < k {
        return Ok(Vec::new());
    }
    Ok((0..=codes.len() - k)
        .map(|i| {
            let part: Vec<u64> = codes[i..i + k].iter().map(|x| *x as u64).collect();
            (i, hash_kmer(&part, 0xABCDEF))
        })
        .collect())
}

pub fn spaced_kmers(seq: &str, pattern: &str, k: usize) -> Result<Vec<(usize, u128)>> {
    if pattern.is_empty() || !pattern.chars().all(|c| matches!(c, '0' | '1')) {
        bail!("spaced pattern 必须是非空 0/1 字符串");
    }
    let ones: Vec<usize> = pattern.match_indices('1').map(|x| x.0).collect();
    if ones.is_empty() {
        bail!("spaced pattern 必须至少包含一个 1");
    }
    if k != ones.len() {
        bail!("k 必须等于 pattern 中 1 的数量");
    }
    let codes = encode(seq);
    if codes.len() < pattern.len() {
        return Ok(Vec::new());
    }
    Ok((0..=codes.len() - pattern.len())
        .map(|i| {
            let part: Vec<u64> = ones.iter().map(|off| codes[i + off] as u64).collect();
            (i, hash_kmer(&part, 0xABCDEF))
        })
        .collect())
}

pub fn minimizers(seq: &str, k: usize, w: usize) -> Result<BTreeSet<u128>> {
    if w == 0 {
        bail!("w 必须为正数");
    }
    let km = kmers(seq, k)?;
    if km.len() < w {
        return Ok(BTreeSet::new());
    }
    Ok((0..=km.len() - w)
        .filter_map(|i| km[i..i + w].iter().min_by_key(|x| x.1).map(|x| x.1))
        .collect())
}

/// 按显式 Open/Closed 定义提取 syncmer；哈希平局采用 first-min。
pub fn syncmers_with_mode(
    seq: &str,
    k: usize,
    s: usize,
    mode: SyncmerMode,
) -> Result<BTreeSet<u128>> {
    if k == 0 || s == 0 || s >= k {
        bail!("要求 0 < s < k");
    }
    if let SyncmerMode::Open { offset } = mode {
        if offset > k - s {
            bail!("open syncmer offset 不能超过 k-s");
        }
    }
    let codes = encode(seq);
    if codes.len() < k {
        return Ok(BTreeSet::new());
    }
    let mut out = BTreeSet::new();
    for i in 0..=codes.len() - k {
        let smers: Vec<u128> = (0..=k - s)
            .map(|j| {
                let part: Vec<u64> = codes[i + j..i + j + s].iter().map(|x| *x as u64).collect();
                hash_kmer(&part, 0xABCDEF)
            })
            .collect();
        let min = *smers.iter().min().unwrap();
        let t = smers.iter().position(|x| *x == min).unwrap();
        let selected = match mode {
            SyncmerMode::Closed => t == 0 || t == k - s,
            SyncmerMode::Open { offset } => t == offset,
        };
        if selected {
            let part: Vec<u64> = codes[i..i + k].iter().map(|x| *x as u64).collect();
            out.insert(hash_kmer(&part, 0xABCDEF));
        }
    }
    Ok(out)
}

/// 兼容历史布尔接口：false 为 Closed，true 为 Open(offset=0)。
#[deprecated(note = "请改用 syncmers_with_mode 和显式 SyncmerMode")]
pub fn syncmers(seq: &str, k: usize, s: usize, closed: bool) -> Result<BTreeSet<u128>> {
    let mode = if closed {
        SyncmerMode::Open { offset: 0 }
    } else {
        SyncmerMode::Closed
    };
    syncmers_with_mode(seq, k, s, mode)
}

pub fn anchor_candidates(seq: &str, k: usize, w1: usize) -> Result<Vec<(usize, u128)>> {
    if w1 == 0 {
        bail!("w1 必须为正数");
    }
    let km = kmers(seq, k)?;
    if km.len() < w1 {
        return Ok(Vec::new());
    }
    let mut seen = BTreeSet::new();
    let mut anchors = Vec::new();
    for i in 0..=km.len() - w1 {
        let best = *km[i..i + w1].iter().min_by_key(|x| x.1).unwrap();
        if seen.insert(best.1) {
            anchors.push(best);
        }
    }
    Ok(anchors)
}

/// 简化 fixed-offset strobemer；不是标准 randstrobe/minstrobe。
pub fn strobemers(seq: &str, k: usize, w1: usize, w2: usize, t: usize) -> Result<BTreeSet<u128>> {
    if t == 0 {
        bail!("t 必须为正数");
    }
    if w1 == 0 || w2 < w1 {
        bail!("要求 0 < w1 <= w2");
    }
    let anchors = anchor_candidates(seq, k, w1)?;
    let mut out = BTreeSet::new();
    for i in 0..anchors.len() {
        let Some(j) = i.checked_add(t).filter(|j| *j < anchors.len()) else {
            continue;
        };
        let (p1, h1) = anchors[i];
        let (p2, h2) = anchors[j];
        let d = p2 - p1;
        if (w1..=w2).contains(&d) {
            out.insert(hash_kmer(&[h1 as u64, h2 as u64], 0x51ED));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classes_and_residues_are_explicit() {
        let map = class_map();
        assert_eq!(
            AA20.chars()
                .map(|c| map[&c.to_string()].clone())
                .collect::<BTreeSet<_>>()
                .len(),
            5
        );
        assert_eq!(clean_protein("ACDU*ac-.", SecPolicy::C), "ACDCXAC");
        assert_eq!(clean_protein("U", SecPolicy::X), "X");
    }

    #[test]
    fn validates_seed_parameters() {
        assert!(kmers("ACD", 0).is_err());
        assert!(spaced_kmers("ACD", "", 0).is_err());
        assert!(spaced_kmers("ACD", "101", 2).is_ok());
        assert!(spaced_kmers("ACD", "101", 3).is_err());
        assert!(spaced_kmers("ACDE", "0000", 0).is_err());
        assert!(syncmers_with_mode("ACD", 3, 3, SyncmerMode::Closed).is_err());
        assert!(strobemers("ACD", 2, 1, 5, 0).is_err());
    }
}
