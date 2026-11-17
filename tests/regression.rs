use chainmer::core::{chainmer_keys, chainmer_keys_legacy, ChainmerParams};
use chainmer::index::{database_keys, SeedIndex};
use chainmer::scoring::weighted_score;
use chainmer::search::{search, search_count, search_with_keys};
use chainmer::seeds::{
    class_map, clean_protein, encode, hash_kmer, kmers, strobemers, syncmers_with_mode,
    to_class_string, SecPolicy, SyncmerMode, AA20,
};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;

const LONG: &str = "MVKAVLAVLAKGGTIAAQLSPEHMPVRMRMLDNNEAAELKALGVPFGRIRLPGGQVAEPIKVDANGVITDAQSSSLIFRGFERILQHQDPRDAAFLTQRICGVCPLSHGLTATNALDELYGVAEHVPKDGYKVTYAPQISDEEVKQMLEDLKRAG";
const FRH: &str = "MTNKKAAAWLVGAIVGGLGMAEAKLDDKAEAKKAAKEAKKAAELKALGVPFGRIRLPGGQVAEPIKVDANGVITDAQSSSLIFRGFERILQHQDPRDAAWQISNEEVKQMLEDLKRAG";

fn coded(text: &str) -> Vec<u64> {
    encode(text).into_iter().map(u64::from).collect()
}

#[test]
fn python_prototype_hash_and_key_parity() {
    let hash = hash_kmer(&coded("ACDEFGH"), 0xABCDEF);
    assert_eq!(hash, 0xd4c2adf8a9645bcc77bdedcf9999e8d2);
    // 历史对照钉住 v0.4 原型参数 delta=15，不随新默认值（delta=5）漂移。
    let legacy_params = ChainmerParams { delta: 15, ..ChainmerParams::default() };
    let keys = chainmer_keys_legacy(LONG, &legacy_params).unwrap();
    assert_eq!(keys.len(), 306);
    assert_eq!(keys.first(), Some(&0x0a1cbb4e9e2bd0adfd911785526697d));
    assert_eq!(keys.last(), Some(&0xffc7a9e54c658efa18c5b069a279c357));
}

#[test]
fn classes_cleaning_and_short_inputs_are_consistent() {
    let map = class_map();
    assert_eq!(AA20.len(), 20);
    assert!(AA20.chars().all(|aa| map.contains_key(&aa.to_string())));
    assert_eq!(map["G"], "b");
    assert_eq!(map["X"], "f");
    assert_eq!(to_class_string("AVLIM"), "aaaaa");
    assert_eq!(clean_protein("AC-D.F*G", SecPolicy::C), "ACDFXG");
    assert_eq!(clean_protein("ACDUFG", SecPolicy::C), "ACDCFG");
    assert_eq!(clean_protein("ACDUFG", SecPolicy::X), "ACDXFG");
    assert!(kmers("", 5).unwrap().is_empty());
    assert!(kmers("ACD", 5).unwrap().is_empty());
    assert!(syncmers_with_mode("ACD", 5, 2, SyncmerMode::Closed)
        .unwrap()
        .is_empty());
    assert!(chainmer_keys("ACD", &ChainmerParams::default())
        .unwrap()
        .is_empty());
}

fn reference_syncmers(seq: &str, k: usize, s: usize, mode: SyncmerMode) -> BTreeSet<u128> {
    let codes = coded(seq);
    let mut expected = BTreeSet::new();
    for i in 0..=codes.len() - k {
        let smers: Vec<u128> = (0..=k - s)
            .map(|j| hash_kmer(&codes[i + j..i + j + s], 0xABCDEF))
            .collect();
        let minimum = *smers.iter().min().unwrap();
        let first_position = smers.iter().position(|value| *value == minimum).unwrap();
        let selected = match mode {
            SyncmerMode::Closed => first_position == 0 || first_position == k - s,
            SyncmerMode::Open { offset } => first_position == offset,
        };
        if selected {
            expected.insert(hash_kmer(&codes[i..i + k], 0xABCDEF));
        }
    }
    expected
}

#[test]
#[allow(deprecated)]
fn syncmer_open_and_closed_match_independent_definitions() {
    for mode in [SyncmerMode::Closed, SyncmerMode::Open { offset: 1 }] {
        assert_eq!(
            syncmers_with_mode(LONG, 5, 2, mode).unwrap(),
            reference_syncmers(LONG, 5, 2, mode)
        );
    }
    let full = |text| hash_kmer(&coded(text), 0xABCDEF);
    // 哈希最小值分别位于内部、左端和右端，不能按氨基酸字典序推断。
    assert!(syncmers_with_mode("AACCC", 5, 2, SyncmerMode::Closed)
        .unwrap()
        .is_empty());
    assert_eq!(
        syncmers_with_mode("AACCC", 5, 2, SyncmerMode::Open { offset: 1 }).unwrap(),
        BTreeSet::from([full("AACCC")])
    );
    assert_eq!(
        syncmers_with_mode("CCCAA", 5, 2, SyncmerMode::Closed).unwrap(),
        BTreeSet::from([full("CCCAA")])
    );
    assert_eq!(
        syncmers_with_mode("CCAAC", 5, 2, SyncmerMode::Closed).unwrap(),
        BTreeSet::from([full("CCAAC")])
    );
    assert!(
        syncmers_with_mode("CCAAC", 5, 2, SyncmerMode::Open { offset: 1 })
            .unwrap()
            .is_empty()
    );
    assert!(
        syncmers_with_mode("CCAAC", 5, 2, SyncmerMode::Open { offset: 2 })
            .unwrap()
            .is_empty()
    );
    // AA 同时位于左右端点且同为最小值；first-min 只把左端视为最小位置。
    assert_eq!(
        syncmers_with_mode("AADAA", 5, 2, SyncmerMode::Open { offset: 0 }).unwrap(),
        BTreeSet::from([full("AADAA")])
    );
    assert!(
        syncmers_with_mode("AADAA", 5, 2, SyncmerMode::Open { offset: 3 })
            .unwrap()
            .is_empty()
    );
    assert!(syncmers_with_mode("ACDE", 5, 2, SyncmerMode::Open { offset: 4 }).is_err());
    // 旧布尔包装保持历史映射，但不再用于 CLI 标签。
    assert_eq!(
        chainmer::seeds::syncmers("CCCAA", 5, 2, false).unwrap(),
        BTreeSet::from([full("CCCAA")])
    );
}

#[test]
fn self_match_cross_family_and_gap_modes() {
    let params = ChainmerParams::default();
    let own = chainmer_keys(LONG, &params).unwrap();
    let other = chainmer_keys(FRH, &params).unwrap();
    assert!(!own.is_empty());
    assert!(own.intersection(&other).count() < own.len());
    let mut no_gap = params.clone();
    no_gap.set_use_gap(false);
    assert_ne!(own, chainmer_keys(LONG, &no_gap).unwrap());
}

#[test]
fn fixed_offset_strobemer_supports_offsets_above_five() {
    let keys = strobemers(LONG, 4, 6, 300, 6).unwrap();
    assert!(!keys.is_empty());
}

#[test]
fn negative_rejection_threshold_is_strictly_greater_than_limit() {
    let params = ChainmerParams::default();
    let query_keys = chainmer_keys(LONG, &params).unwrap();
    let negative: BTreeSet<u128> = query_keys.iter().copied().take(3).collect();
    let db = BTreeMap::from([
        (
            "A".to_string(),
            query_keys.difference(&negative).copied().collect(),
        ),
        ("B".to_string(), query_keys.clone()),
    ]);
    for (limit, contains_b) in [(0, false), (2, false), (3, true)] {
        let hits: BTreeMap<String, usize> =
            search_with_keys(LONG, &db, &params, 1, &negative, limit)
                .unwrap()
                .into_iter()
                .collect();
        assert!(hits.contains_key("A"));
        assert_eq!(hits.contains_key("B"), contains_b);
    }
}

#[test]
fn weighted_score_orders_positive_enriched_candidate_first() {
    let weights = chainmer::contrastive_weights(&[BTreeSet::from([1, 2])], &[BTreeSet::from([3])]);
    let query = BTreeSet::from([1, 9]);
    let positive = weighted_score(&query, &BTreeSet::from([1]), &weights, 0.0).unwrap();
    let negative = weighted_score(&query, &BTreeSet::from([3]), &weights, 0.0).unwrap();
    assert!(positive > negative);
}

#[test]
fn index_roundtrip_preserves_metadata_postings_and_results() {
    let db = BTreeMap::from([
        ("a".to_string(), LONG.to_string()),
        ("b".to_string(), LONG.to_string()),
    ]);
    let params = ChainmerParams::default();
    let index = SeedIndex::build(&db, params.clone()).unwrap();
    let path =
        std::env::temp_dir().join(format!("chainmer-regression-{}.json", std::process::id()));
    index.save(&path).unwrap();
    let loaded = SeedIndex::load(&path).unwrap();
    fs::remove_file(&path).unwrap();
    assert_eq!(loaded, index);
    assert_eq!(loaded.meta.db.n_sequences, 2);
    assert_eq!(loaded.meta.class_map["G"], "b");
    assert!(loaded.meta.hash_scheme.contains("64-effective"));
    let keys = database_keys(&db, &params).unwrap();
    let indexed = search_count(LONG, &loaded, &params, 1).unwrap();
    let brute = search_with_keys(LONG, &keys, &params, 1, &BTreeSet::new(), 0).unwrap();
    assert_eq!(indexed, brute);
}

#[test]
fn chain_score_penalizes_extra_unrelated_target_keys() {
    let elongated = format!("{LONG}{FRH}");
    let db = BTreeMap::from([
        ("self".to_string(), LONG.to_string()),
        ("elongated".to_string(), elongated),
        ("other".to_string(), FRH.to_string()),
    ]);
    let params = ChainmerParams::default();
    let index = SeedIndex::build(&db, params.clone()).unwrap();
    let hits = search(LONG, &index, &params, 1, &BTreeSet::new(), 0, None).unwrap();
    let own = hits.iter().find(|hit| hit.id == "self").unwrap();
    let long = hits.iter().find(|hit| hit.id == "elongated").unwrap();
    assert_eq!(own.shared, long.shared);
    assert!(own.weighted_jaccard > long.weighted_jaccard);
    assert!(own.score > long.score);
}

#[test]
fn database_frequency_downweights_common_keys() {
    let db = BTreeMap::from([
        ("a".to_string(), LONG.to_string()),
        ("b".to_string(), FRH.to_string()),
    ]);
    let params = ChainmerParams::default();
    let index = SeedIndex::build(&db, params.clone()).unwrap();
    // IDF 只定义在索引覆盖的锚点对键上（密集间距键不进入 CSR）
    let a = chainmer::core::chainmer_pair_keys(LONG, &params).unwrap();
    let b = chainmer::core::chainmer_pair_keys(FRH, &params).unwrap();
    let common = a.intersection(&b).next().unwrap();
    let rare = a.difference(&b).next().unwrap();
    assert!(index.idf(common) < index.idf(rare));
}

#[test]
fn adaptive_mode_does_not_change_native_seed_keys() {
    // adaptive 仅新增排序层；原生链式键必须与历史行为逐项一致
    let params = ChainmerParams::default();
    let keys = chainmer_keys(LONG, &params).unwrap();
    assert!(!keys.is_empty());
    let legacy = chainmer_keys_legacy(LONG, &ChainmerParams { delta: 15, ..params }).unwrap();
    assert_eq!(legacy.len(), 306);
    // 子序列通道与链式键空间不相交（blake2b-64 低 64 位 vs 双 splitmix u128）
    let sub = chainmer::adaptive::subset_keys(LONG);
    assert!(!sub.is_empty());
    let sub_as_u128: BTreeSet<u128> = sub.iter().copied().collect();
    assert!(keys.is_disjoint(&sub_as_u128));
}

#[test]
fn adaptive_family_context_and_decoupled_gate() {
    use chainmer::adaptive::{adaptive_search, AdaptiveParams, FamilyContext, Regime};
    let mut para = LONG.to_string();
    para.replace_range(40..60, "AAAAAAAAAAAAAAAAAAAA");
    let ids = vec![
        "SpA__q".to_string(),
        "SpA__para".to_string(),
        "SpB__orth".to_string(),
        "SpC__orth".to_string(),
    ];
    let seqs = vec![LONG.to_string(), para, LONG.to_string(), FRH.to_string()];
    let cp = ChainmerParams::default();
    let ap = AdaptiveParams::default();
    let ctx = FamilyContext::build(ids.clone(), &seqs, &ap, &cp).unwrap();
    let r = adaptive_search(&ctx, "SpA__q", LONG, &ap, &cp).unwrap();
    assert_eq!(r.hits.len(), 3);
    // 召回门始终为链式 Jaccard：自身物种旁系与远缘直系都不得被背景扣除否决
    assert!(r.hits.iter().all(|h| h.gate >= 0.0));
    let gate_of = |id: &str| r.hits.iter().find(|h| h.id == id).unwrap().gate;
    assert!(gate_of("SpB__orth") > 0.0);
    // 体制标签合法且与家族门控一致
    if ctx.fam_divergent {
        assert_eq!(r.regime, Regime::Divergent);
        assert!(r.hits.iter().all(|h| (h.score - h.gate).abs() < 1e-12));
    }
}

#[test]
fn index_v7_adaptive_layer_roundtrip() {
    use chainmer::adaptive::{family_divergent, spacing_hist, subset_keys};
    use chainmer::core::{chainmer_dense_keys, chainmer_pair_keys};
    let mut para = LONG.to_string();
    para.replace_range(40..60, "AAAAAAAAAAAAAAAAAAAA");
    let db = BTreeMap::from([
        ("SpA__q".to_string(), LONG.to_string()),
        ("SpA__para".to_string(), para),
        ("SpB__orth".to_string(), FRH.to_string()),
    ]);
    let params = ChainmerParams::default();
    let index = SeedIndex::build(&db, params.clone()).unwrap();
    // 判别层与独立重算逐项一致（链键表按锚点对键、dense 表按间距键）
    let hists: Vec<Vec<f64>> = db.values().map(|s| spacing_hist(s)).collect();
    for (tid, (id, seq)) in db.iter().enumerate() {
        let ck: BTreeSet<u128> = index.chain_keys[tid].iter().copied().collect();
        assert_eq!(ck, chainmer_pair_keys(seq, &params).unwrap(), "{id}");
        let dk: BTreeSet<u128> = index.dense_keys[tid].iter().copied().collect();
        assert_eq!(dk, chainmer_dense_keys(seq, &params).unwrap(), "{id}");
        let sk: BTreeSet<u128> = index.subset_keys[tid].iter().map(|&k| k as u128).collect();
        assert_eq!(sk, subset_keys(seq));
        assert_eq!(index.hists[tid], hists[tid]);
    }
    assert_eq!(
        index.meta.fam_divergent,
        family_divergent(&hists, index.meta.div_threshold)
    );
    // 二进制 roundtrip 保留全部判别层字段
    let path = std::env::temp_dir().join(format!(
        "chainmer-regression-v7-{}.ochidx",
        std::process::id()
    ));
    index.save(&path).unwrap();
    let loaded = SeedIndex::load(&path).unwrap();
    fs::remove_file(&path).unwrap();
    assert_eq!(loaded, index);
}

#[test]
fn search_adaptive_index_path_matches_brute_force() {
    use chainmer::adaptive::{adaptive_search, AdaptiveParams, FamilyContext};
    use chainmer::search::search_adaptive;
    let mut para = LONG.to_string();
    para.replace_range(40..60, "AAAAAAAAAAAAAAAAAAAA");
    let ids = vec![
        "SpA__q".to_string(),
        "SpA__para".to_string(),
        "SpB__orth".to_string(),
        "SpC__orth".to_string(),
    ];
    let seqs = vec![LONG.to_string(), para, LONG.to_string(), FRH.to_string()];
    let db: BTreeMap<String, String> = ids.iter().cloned().zip(seqs.iter().cloned()).collect();
    let cp = ChainmerParams::default();
    let ap = AdaptiveParams::default();
    let ctx = FamilyContext::build(ids.clone(), &seqs, &ap, &cp).unwrap();
    let index = SeedIndex::build(&db, cp.clone()).unwrap();
    for qid in &ids {
        let brute = adaptive_search(&ctx, qid, &db[qid], &ap, &cp).unwrap();
        let indexed = search_adaptive(&db[qid], qid, &index, &cp, &ap, 1).unwrap();
        assert_eq!(brute.regime, indexed.regime);
        assert_eq!(brute.fam_divergent, indexed.fam_divergent);
        assert_eq!(brute.hits.len(), indexed.hits.len());
        // 逐 ID 得分/召回门完全一致（平局裁决序不同：文件序 vs catalog 序）
        for hit in &brute.hits {
            let other = indexed.hits.iter().find(|h| h.id == hit.id).unwrap();
            assert!((hit.score - other.score).abs() < 1e-12);
            assert!((hit.gate - other.gate).abs() < 1e-12);
        }
        // 无平局时排序也完全一致
        let mut scores: Vec<f64> = brute.hits.iter().map(|h| h.score).collect();
        scores.sort_by(|a, b| b.total_cmp(a));
        let has_tie = scores.windows(2).any(|w| w[0] == w[1]);
        if !has_tie {
            let brute_order: Vec<&str> = brute.hits.iter().map(|h| h.id.as_str()).collect();
            let indexed_order: Vec<&str> = indexed.hits.iter().map(|h| h.id.as_str()).collect();
            assert_eq!(brute_order, indexed_order);
        }
    }
}
