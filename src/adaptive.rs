//! 体制自适应判别模式（adaptive，V10 方案的 Rust 实现）。
//!
//! 链式种子召回层完全不变：候选接受门仍是共线键 IDF-Jaccard `Jc > 0`。
//! 本模块只在保守体制下对已召回候选做判别式重排：
//!   score = Jc + β·Js − λ·Jbg
//! 其中 Js 为诊断子序列通道（span20/weight8 四掩码、步长 2），Jbg 为判别式
//! 背景扣除（背景＝每多拷贝物种中较不相似于查询的拷贝，无监督旁系代理）。
//! 体制门控：同源集 RAW 锚点间距直方图的最小两两余弦 < div_threshold ⇒ 分歧
//! 体制，走原生打分；保守体制再按旁系相似度比（par_j/top1_j > 0.7）决定是否
//! 叠加背景扣除。全部信号仅来自库内序列，不需要树/HMM/注释等额外输入。
//!
//! 口径与 Python 原型 proto_adapt.py 逐项对齐：blake2b-64 子序列哈希、
//! IDF=log1p(n/df) 按候选库（排除查询自身）统计、稳定排序。
//! 评分采用 f64 精度：旧基准管线 scipy float32 CSR 在小分值下有显著舍入
//! 伪影（相对误差可达数十％），adaptive 不复制该伪影，对照实验须同用 f64。

use crate::core::{chainmer_keys, ChainmerParams};
use crate::seeds::{clean_protein, SecPolicy};
use anyhow::Result;
use blake2::digest::{consts::U8, Digest};
use blake2::Blake2b;
use std::collections::{BTreeMap, BTreeSet};

type Blake2b64 = Blake2b<U8>;

/// 背景扣除强度 λ（原型实测默认）。
pub const DEFAULT_LAMBDA: f64 = 1.0;
/// 子序列通道权重 β（原型实测默认）。
pub const DEFAULT_BETA: f64 = 1.0;
/// 家族级架构门控阈值：最小两两余弦低于此值判分歧体制。
pub const DEFAULT_DIV_THRESHOLD: f64 = 0.8;
/// 高旁系风险阈值：旁系代理/top1 相似度比高于此值才叠加背景扣除。
pub const RISK_RATIO: f64 = 0.7;

/// 非连续 spaced 掩码（span=20, weight=8），采样散在枝诊断位点。
pub const MASKS: [[usize; 8]; 4] = [
    [0, 2, 4, 6, 9, 12, 15, 19],
    [0, 3, 5, 8, 10, 13, 16, 19],
    [1, 3, 6, 8, 11, 14, 17, 19],
    [0, 2, 5, 7, 10, 12, 15, 18],
];
const SPAN: usize = 20;
const STRIDE: usize = 2;

/// 与原型 hkey 一致：blake2b-64，残基字符串大端解释为 u64。
fn hkey(chars: &[u8]) -> u64 {
    let mut hasher = Blake2b64::new();
    hasher.update(chars);
    let digest = hasher.finalize();
    u64::from_be_bytes(digest.into())
}

fn clean(seq: &str) -> String {
    clean_protein(seq, SecPolicy::C)
}

/// 诊断子序列键集：span20 窗口、步长 2、四掩码；覆盖 X 的采样跳过。
pub fn subset_keys(seq: &str) -> BTreeSet<u128> {
    let s = clean(seq);
    let b = s.as_bytes();
    let mut out = BTreeSet::new();
    if b.len() < SPAN {
        return out;
    }
    let mut i = 0;
    while i + SPAN <= b.len() {
        let win = &b[i..i + SPAN];
        for mask in &MASKS {
            let t: Vec<u8> = mask.iter().map(|p| win[*p]).collect();
            if t.contains(&b'X') {
                continue;
            }
            out.insert(hkey(&t) as u128);
        }
        i += STRIDE;
    }
    out
}

/// 升序去重 Vec 版（免 BTreeSet 节点分配；检索热路径用）。
pub fn subset_keys_vec(seq: &str) -> Vec<u128> {
    let s = clean(seq);
    let b = s.as_bytes();
    let mut out: Vec<u128> = Vec::new();
    if b.len() < SPAN {
        return out;
    }
    let mut i = 0;
    while i + SPAN <= b.len() {
        let win = &b[i..i + SPAN];
        for mask in &MASKS {
            let t: Vec<u8> = mask.iter().map(|p| win[*p]).collect();
            if t.contains(&b'X') {
                continue;
            }
            out.push(hkey(&t) as u128);
        }
        i += STRIDE;
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// RAW 通道 minimizer 锚点→锚点对间距桶直方图（架构签名，归一化）。
///
/// 与原型 spacing_hist 逐项一致：k=5、w=6、d∈[30,300]、δ=5、每锚点向后
/// 至多 12 对；含 X 的 k-mer 不参与窗口最小值；相邻重复锚点位置去重。
pub fn spacing_hist(seq: &str) -> Vec<f64> {
    const K: usize = 5;
    const W: usize = 6;
    const DMIN: usize = 30;
    const DMAX: usize = 300;
    const DELTA: usize = 5;
    let nb = (DMAX - DMIN) / DELTA + 1;
    let s = clean(seq);
    let b = s.as_bytes();
    let mut hist = vec![0.0f64; nb];
    if b.len() < K {
        return hist;
    }
    let hashes: Vec<Option<u64>> = (0..=b.len() - K)
        .map(|i| {
            let part = &b[i..i + K];
            if part.contains(&b'X') {
                None
            } else {
                Some(hkey(part))
            }
        })
        .collect();
    let mut anchors: Vec<usize> = Vec::new();
    for i in 0..=hashes.len().saturating_sub(W) {
        let mut best: Option<(u64, usize)> = None;
        for (j, h) in hashes[i..i + W].iter().enumerate() {
            if let Some(h) = h {
                if best.map(|(bh, _)| *h < bh).unwrap_or(true) {
                    best = Some((*h, i + j));
                }
            }
        }
        if let Some((_, p)) = best {
            if anchors.last() != Some(&p) {
                anchors.push(p);
            }
        }
    }
    for a in 0..anchors.len() {
        for b_idx in (a + 1)..(a + 13).min(anchors.len()) {
            let d = anchors[b_idx] - anchors[a];
            if (DMIN..=DMAX).contains(&d) {
                hist[(d - DMIN) / DELTA] += 1.0;
            }
        }
    }
    let tot: f64 = hist.iter().sum();
    if tot > 0.0 {
        for v in &mut hist {
            *v /= tot;
        }
    }
    hist
}

pub fn cosine(a: &[f64], b: &[f64]) -> f64 {
    let na: f64 = a.iter().map(|x| x * x).sum::<f64>().sqrt();
    let nb: f64 = b.iter().map(|x| x * x).sum::<f64>().sqrt();
    if na > 0.0 && nb > 0.0 {
        a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>() / (na * nb)
    } else {
        0.0
    }
}

/// 家族级架构门控：同源集间距直方图最小两两余弦 < threshold ⇒ 分歧体制。
pub fn family_divergent(hists: &[Vec<f64>], threshold: f64) -> bool {
    for i in 0..hists.len() {
        for j in (i + 1)..hists.len() {
            if cosine(&hists[i], &hists[j]) < threshold {
                return true;
            }
        }
    }
    false
}

/// IDF 加权 Jaccard，与 bench.similarity_matrix df_from="db" 同定义：
/// df 按评分时的候选库统计；全链路 f64 精度（不复制旧管线 float32 舍入伪影）。
pub struct IdfJaccard {
    idf: BTreeMap<u128, f64>,
    db_weights: Vec<f64>,
}

impl IdfJaccard {
    pub fn new(db: &[BTreeSet<u128>]) -> Self {
        let mut df: BTreeMap<u128, usize> = BTreeMap::new();
        for keys in db {
            for k in keys {
                *df.entry(*k).or_default() += 1;
            }
        }
        let n = db.len().max(1) as f64;
        let idf: BTreeMap<u128, f64> = df
            .into_iter()
            .map(|(k, d)| (k, (1.0 + n / d.max(1) as f64).ln()))
            .collect();
        let db_weights = db.iter().map(|keys| Self::weight_sum(keys, &idf)).collect();
        Self { idf, db_weights }
    }

    fn weight_sum(keys: &BTreeSet<u128>, idf: &BTreeMap<u128, f64>) -> f64 {
        keys.iter().filter_map(|k| idf.get(k)).sum()
    }

    /// 查询键集对每个库序列的相似度。
    pub fn scores(&self, query: &BTreeSet<u128>, db: &[BTreeSet<u128>]) -> Vec<f64> {
        let qw = Self::weight_sum(query, &self.idf);
        db.iter()
            .zip(&self.db_weights)
            .map(|(d, dw)| {
                let inter: f64 = query
                    .intersection(d)
                    .filter_map(|k| self.idf.get(k))
                    .sum();
                let denom = qw + dw - inter;
                if denom > 0.0 {
                    (inter / denom).clamp(0.0, 1.0)
                } else {
                    0.0
                }
            })
            .collect()
    }
}

/// 序列 ID 的物种前缀（species__accession 约定）。
pub fn species_of(id: &str) -> &str {
    id.split("__").next().unwrap_or(id)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Regime {
    /// 分歧体制：原生 Chainmer 打分（保护架构分歧甜区）。
    Divergent,
    /// 保守＋高旁系风险：Jc + β·Js − λ·Jbg。
    ConservedBg,
    /// 保守＋低风险：Jc + β·Js（保护饱和家族）。
    ConservedSub,
}

impl Regime {
    pub fn label(self) -> &'static str {
        match self {
            Self::Divergent => "divergent",
            Self::ConservedBg => "conserved-bg",
            Self::ConservedSub => "conserved-sub",
        }
    }
}

pub struct AdaptiveParams {
    pub lambda: f64,
    pub beta: f64,
    pub div_threshold: f64,
}

impl Default for AdaptiveParams {
    fn default() -> Self {
        Self {
            lambda: DEFAULT_LAMBDA,
            beta: DEFAULT_BETA,
            div_threshold: DEFAULT_DIV_THRESHOLD,
        }
    }
}

/// 家族级预计算：全同源集（含查询）的键、子序列键、架构直方图与体制门控。
/// ID 顺序必须为 FASTA 文件顺序（平局裁决与原型一致）。
pub struct FamilyContext {
    pub ids: Vec<String>,
    pub ck: Vec<BTreeSet<u128>>,
    pub sk: Vec<BTreeSet<u128>>,
    pub hists: Vec<Vec<f64>>,
    pub fam_divergent: bool,
}

impl FamilyContext {
    pub fn build(
        ids: Vec<String>,
        seqs: &[String],
        ap: &AdaptiveParams,
        cp: &ChainmerParams,
    ) -> Result<Self> {
        let ck: Vec<BTreeSet<u128>> = seqs
            .iter()
            .map(|s| chainmer_keys(s, cp))
            .collect::<Result<_>>()?;
        let sk: Vec<BTreeSet<u128>> = seqs.iter().map(|s| subset_keys(s)).collect();
        let hists: Vec<Vec<f64>> = seqs.iter().map(|s| spacing_hist(s)).collect();
        let fam_divergent = family_divergent(&hists, ap.div_threshold);
        Ok(Self { ids, ck, sk, hists, fam_divergent })
    }

    fn index_of(&self, id: &str) -> Option<usize> {
        self.ids.iter().position(|x| x == id)
    }
}

pub struct AdaptiveHit {
    pub id: String,
    /// 判别排序得分（体制门控后）。
    pub score: f64,
    /// 召回/接受门：共线键 IDF-Jaccard Jc（与原生一致，>0 即命中）。
    pub gate: f64,
}

pub struct AdaptiveQueryResult {
    pub query: String,
    pub regime: Regime,
    pub fam_divergent: bool,
    pub conserved: bool,
    pub hits: Vec<AdaptiveHit>,
}

/// 对单个查询执行 adaptive 检索；候选库＝全同源集排除查询自身。
///
/// IDF 按候选库（排除查询）统计，与原型 similarity_matrix(df_from="db") 一致。
pub fn adaptive_search(
    ctx: &FamilyContext,
    qid: &str,
    qseq: &str,
    ap: &AdaptiveParams,
    cp: &ChainmerParams,
) -> Result<AdaptiveQueryResult> {
    let (qk, qsk) = match ctx.index_of(qid) {
        Some(i) => (ctx.ck[i].clone(), ctx.sk[i].clone()),
        None => (chainmer_keys(qseq, cp)?, subset_keys(qseq)),
    };
    // 候选库排除查询自身；IDF 词表仍按全集统计
    let db_idx: Vec<usize> = (0..ctx.ids.len()).filter(|&i| ctx.ids[i] != qid).collect();
    let db_ck: Vec<BTreeSet<u128>> = db_idx.iter().map(|&i| ctx.ck[i].clone()).collect();
    let db_sk: Vec<BTreeSet<u128>> = db_idx.iter().map(|&i| ctx.sk[i].clone()).collect();
    let db_ids: Vec<String> = db_idx.iter().map(|&i| ctx.ids[i].clone()).collect();
    Ok(adaptive_core(
        qid,
        &qk,
        &qsk,
        &db_ids,
        &db_ck,
        &db_sk,
        ctx.fam_divergent,
        ap,
    ))
}

/// 判别核心：给定查询键集与候选库视图（已排除查询自身），执行体制门控打分。
///
/// 暴力路径（FamilyContext）与索引路径（search_adaptive）共用此函数，保证
/// 两条路径数值完全一致。db_ids 的顺序即平局裁决顺序（稳定排序）。
pub fn adaptive_core(
    qid: &str,
    qk: &BTreeSet<u128>,
    qsk: &BTreeSet<u128>,
    db_ids: &[String],
    db_ck: &[BTreeSet<u128>],
    db_sk: &[BTreeSet<u128>],
    fam_divergent: bool,
    ap: &AdaptiveParams,
) -> AdaptiveQueryResult {
    let jc_sc = IdfJaccard::new(db_ck);
    let jc = jc_sc.scores(qk, db_ck);
    let js_sc = IdfJaccard::new(db_sk);
    let js = js_sc.scores(qsk, db_sk);

    // 物种多拷贝分组（按首次出现顺序，与原型 dict 插入序一致）
    let mut sp_order: Vec<String> = Vec::new();
    let mut by_sp: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (i, id) in db_ids.iter().enumerate() {
        let sp = species_of(id).to_string();
        if !by_sp.contains_key(&sp) {
            sp_order.push(sp.clone());
        }
        by_sp.entry(sp).or_default().push(i);
    }
    let qsp = species_of(qid);
    let top1_j = jc.iter().cloned().fold(0.0f64, f64::max);

    // 背景集与旁系风险：同物种取其余全部拷贝；他物种取较不相似拷贝（代理）。
    // 多拷贝判定按全同源集（含查询）：查询自身物种只需库内≥1 个其余拷贝。
    let mut bg: BTreeSet<u128> = BTreeSet::new();
    let mut par_j = 0.0f64;
    for sp in &sp_order {
        let cand = &by_sp[sp];
        let is_qsp = sp == qsp;
        if cand.is_empty() || (!is_qsp && cand.len() < 2) {
            continue;
        }
        if is_qsp {
            for &i in cand {
                bg.extend(&db_ck[i]);
                par_j = par_j.max(jc[i]);
            }
        } else {
            // 较不相似拷贝（代理）：jc 平局时按 ID 确定性裁决，保证索引/暴力两路径一致。
            let least = *cand
                .iter()
                .min_by(|&&a, &&b| {
                    jc[a]
                        .total_cmp(&jc[b])
                        .then_with(|| db_ids[a].cmp(&db_ids[b]))
                })
                .unwrap();
            bg.extend(&db_ck[least]);
            par_j = par_j.max(jc[least]);
        }
    }
    let jbg = if bg.is_empty() {
        vec![0.0; db_ids.len()]
    } else {
        jc_sc.scores(&bg, &db_ck)
    };
    let conserved = top1_j > 0.0 && par_j / top1_j > RISK_RATIO;

    let regime = if fam_divergent {
        Regime::Divergent
    } else if conserved {
        Regime::ConservedBg
    } else {
        Regime::ConservedSub
    };
    let scores: Vec<f64> = match regime {
        Regime::Divergent => jc.clone(),
        Regime::ConservedBg => jc
            .iter()
            .zip(&js)
            .zip(&jbg)
            .map(|((c, s), b)| c + ap.beta * s - ap.lambda * b)
            .collect(),
        Regime::ConservedSub => jc.iter().zip(&js).map(|(c, s)| c + ap.beta * s).collect(),
    };

    // 稳定排序：同分保持传入的 db_ids 顺序（暴力路径为 FASTA 序，索引路径为 catalog 序）
    let mut order: Vec<usize> = (0..db_ids.len()).collect();
    order.sort_by(|&a, &b| scores[b].total_cmp(&scores[a]));
    let hits = order
        .into_iter()
        .map(|i| AdaptiveHit {
            id: db_ids[i].clone(),
            score: scores[i],
            gate: jc[i],
        })
        .collect();
    AdaptiveQueryResult {
        query: qid.to_string(),
        regime,
        fam_divergent,
        conserved,
        hits,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const LONG: &str = "MVKAVLAVLAKGGTIAAQLSPEHMPVRMRMLDNNEAAELKALGVPFGRIRLPGGQVAEPIKVDANGVITDAQSSSLIFRGFERILQHQDPRDAAFLTQRICGVCPLSHGLTATNALDELYGVAEHVPKDGYKVTYAPQISDEEVKQMLEDLKRAG";

    #[test]
    fn hkey_matches_blake2b64_big_endian() {
        let mut hasher = Blake2b64::new();
        hasher.update(b"ACDEFGHI");
        let expect = u64::from_be_bytes(hasher.finalize().into());
        assert_eq!(hkey(b"ACDEFGHI"), expect);
    }

    #[test]
    fn subset_keys_mask_x_and_require_span() {
        assert!(subset_keys("ACDEFG").is_empty());
        let keys = subset_keys(LONG);
        assert!(!keys.is_empty());
        assert!(subset_keys(&"X".repeat(60)).is_empty());
        // 掩码覆盖 X 的采样被跳过：键数不多于无 X 对照
        let with_x = format!("{}X{}", &LONG[..80], &LONG[81..]);
        assert!(subset_keys(&with_x).len() <= keys.len());
    }

    #[test]
    fn spacing_hist_normalized_and_empty_safe() {
        let h = spacing_hist(LONG);
        assert_eq!(h.len(), 55);
        let tot: f64 = h.iter().sum();
        assert!((tot - 1.0).abs() < 1e-9 || tot == 0.0);
        assert_eq!(spacing_hist("ACD").iter().sum::<f64>(), 0.0);
        assert!((cosine(&h, &h) - 1.0).abs() < 1e-12);
        assert_eq!(cosine(&vec![0.0; 55], &h), 0.0);
    }

    #[test]
    fn family_gate_separates_clusters() {
        let a = vec![1.0, 0.0, 0.0];
        let b = vec![0.9, 0.1, 0.0];
        let c = vec![0.0, 0.0, 1.0];
        assert!(!family_divergent(&[a.clone(), b], 0.8));
        assert!(family_divergent(&[a, c], 0.8));
        assert!(!family_divergent(&[vec![1.0, 0.0]], 0.8));
    }

    #[test]
    fn idf_jaccard_self_one_and_zero_disjoint() {
        let db = vec![
            BTreeSet::from([1u128, 2, 3]),
            BTreeSet::from([2u128, 3, 4]),
        ];
        let sc = IdfJaccard::new(&db);
        let s = sc.scores(&BTreeSet::from([1u128, 2, 3]), &db);
        assert!((s[0] - 1.0).abs() < 1e-12);
        assert!(s[1] > 0.0 && s[1] < 1.0);
        let z = sc.scores(&BTreeSet::from([99u128]), &db);
        assert_eq!(z, vec![0.0, 0.0]);
        // 键 1 仅 1 条 ⇒ idf=ln(3)；f64 精度口径
        let expect = (3.0f64).ln();
        assert!((sc.idf[&1] - expect).abs() < 1e-12);
    }

    #[test]
    fn species_prefix_convention() {
        assert_eq!(
            species_of("Methanosarcina_acetivorans_C2A__AAM07885.1"),
            "Methanosarcina_acetivorans_C2A"
        );
        assert_eq!(species_of("no_separator"), "no_separator");
    }

    #[test]
    fn adaptive_gate_is_always_chainmer_jaccard() {
        // 合成双物种多拷贝库：gate 必须等于 Jc（解耦召回门），排序与门解耦
        let mut s2 = LONG.to_string();
        s2.replace_range(40..60, "AAAAAAAAAAAAAAAAAAAA");
        let ids = vec![
            "SpA__q".to_string(),
            "SpA__para".to_string(),
            "SpB__orth".to_string(),
        ];
        let seqs = vec![LONG.to_string(), s2, LONG.to_string()];
        let cp = ChainmerParams::default();
        let ap = AdaptiveParams::default();
        let ctx = FamilyContext::build(ids.clone(), &seqs, &ap, &cp).unwrap();
        let r = adaptive_search(&ctx, "SpA__q", LONG, &ap, &cp).unwrap();
        let db_ck: Vec<BTreeSet<u128>> = ctx.ck[1..].to_vec();
        let jc = IdfJaccard::new(&db_ck).scores(&ctx.ck[0], &db_ck);
        for hit in &r.hits {
            let i = ids[1..].iter().position(|x| x == &hit.id).unwrap();
            assert!((hit.gate - jc[i]).abs() < 1e-12);
        }
        if r.regime != Regime::Divergent {
            assert!(r.hits.iter().any(|h| (h.score - h.gate).abs() > 1e-9));
        }
    }

    #[test]
    fn divergent_family_keeps_native_scores() {
        // 构造架构分歧家族：一条序列的间距签名与其余显著分离
        let ids = vec!["SpA__q".to_string(), "SpB__x".to_string()];
        let seqs = vec![LONG.to_string(), "A".repeat(200)];
        let cp = ChainmerParams::default();
        let ap = AdaptiveParams::default();
        let ctx = FamilyContext::build(ids, &seqs, &ap, &cp).unwrap();
        let r = adaptive_search(&ctx, "SpA__q", LONG, &ap, &cp).unwrap();
        if ctx.fam_divergent {
            assert_eq!(r.regime, Regime::Divergent);
            for hit in &r.hits {
                assert!((hit.score - hit.gate).abs() < 1e-12);
            }
        }
    }
}