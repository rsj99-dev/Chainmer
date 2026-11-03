use orthochain::core::{orthochain_keys, orthochain_keys_legacy, OrthoParams};
use orthochain::index::{database_keys, SeedIndex};
use orthochain::scoring::weighted_score;
use orthochain::search::{search, search_count, search_with_keys};
use orthochain::seeds::{
    class_map, clean_protein, encode, hash_kmer, kmers, strobemers, syncmers, to_class_string,
    SecPolicy, AA20,
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
    let keys = orthochain_keys_legacy(LONG, &OrthoParams::default()).unwrap();
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
    assert!(syncmers("ACD", 5, 2, false).unwrap().is_empty());
    assert!(orthochain_keys("ACD", &OrthoParams::default())
        .unwrap()
        .is_empty());
}

fn reference_syncmers(seq: &str, k: usize, s: usize, closed: bool) -> BTreeSet<u128> {
    let codes = coded(seq);
    let mut expected = BTreeSet::new();
    for i in 0..=codes.len() - k {
        let smers: Vec<u128> = (0..=k - s)
            .map(|j| hash_kmer(&codes[i + j..i + j + s], 0xABCDEF))
            .collect();
        let minimum = *smers.iter().min().unwrap();
        let position = smers.iter().position(|value| *value == minimum).unwrap();
        if position == 0 || (!closed && position == k - s) {
            expected.insert(hash_kmer(&codes[i..i + k], 0xABCDEF));
        }
    }
    expected
}

#[test]
fn syncmer_open_and_closed_match_definitions() {
    assert_eq!(
        syncmers(LONG, 5, 2, false).unwrap(),
        reference_syncmers(LONG, 5, 2, false)
    );
    assert_eq!(
        syncmers(LONG, 5, 2, true).unwrap(),
        reference_syncmers(LONG, 5, 2, true)
    );
}

#[test]
fn self_match_cross_family_and_gap_modes() {
    let params = OrthoParams::default();
    let own = orthochain_keys(LONG, &params).unwrap();
    let other = orthochain_keys(FRH, &params).unwrap();
    assert!(!own.is_empty());
    assert!(own.intersection(&other).count() < own.len());
    let mut no_gap = params.clone();
    no_gap.set_use_gap(false);
    assert_ne!(own, orthochain_keys(LONG, &no_gap).unwrap());
}

#[test]
fn fixed_offset_strobemer_supports_offsets_above_five() {
    let keys = strobemers(LONG, 4, 6, 300, 6).unwrap();
    assert!(!keys.is_empty());
}

#[test]
fn negative_rejection_threshold_is_strictly_greater_than_limit() {
    let params = OrthoParams::default();
    let query_keys = orthochain_keys(LONG, &params).unwrap();
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
    let weights =
        orthochain::contrastive_weights(&[BTreeSet::from([1, 2])], &[BTreeSet::from([3])]);
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
    let params = OrthoParams::default();
    let index = SeedIndex::build(&db, params.clone()).unwrap();
    let path =
        std::env::temp_dir().join(format!("orthochain-regression-{}.json", std::process::id()));
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
    let params = OrthoParams::default();
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
    let params = OrthoParams::default();
    let index = SeedIndex::build(&db, params.clone()).unwrap();
    let a = orthochain_keys(LONG, &params).unwrap();
    let b = orthochain_keys(FRH, &params).unwrap();
    let common = a.intersection(&b).next().unwrap();
    let rare = a.difference(&b).next().unwrap();
    assert!(index.idf(common) < index.idf(rare));
}
