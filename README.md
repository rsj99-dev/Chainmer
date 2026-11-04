<p align="center">
  <img src="Chainmer_Logo.png" alt="Chainmer Logo" width="420">
</p>

# Chainmer v1.0.0

**Language:** English | [中文](README_CN.md)

Chainmer is a Rust library and command-line tool for **protein homology candidate retrieval and prefiltering**. Version 1.0.0 replaces the former single-track class-equivalence key count with **raw/class dual-track positional seeds, database-calibrated IDF, and collinear chain scoring**. The legacy count and SWPN interfaces remain available as ablation controls.

> A Seed Hit only indicates that candidates share a seed; it does not establish homology or orthology. Orthology requires downstream analyses such as gene tree–species tree reconciliation. HydDB subgroup labels are not equivalent to orthologous groups.

## Core components

1. **Dual-track recall**: the raw track uses original amino acids with `kc=5`, while the class track uses substitution-equivalence classes with `class_kc=8`; the two key types use different hash salts.
2. **Positioned anchor pairs**: both anchor positions are stored for queries and targets, rather than treating seeds as unordered key sets.
3. **Database frequency calibration**: IDF is calculated from document frequencies in the current database, automatically down-weighting frequent keys.
4. **Duplicate contribution constraint**: retained positional hits for the same key within one candidate share a single IDF weight. A stable-hash bottom-k retains at most 64 representative hits per key, making the retained set independent of posting order; hit counts before and after truncation are reported.
5. **Collinear chain scoring**: optimal chains are computed within staggered diagonal bins, with compatibility requiring all four anchor endpoints to increase strictly in both query and target. The final score is the sum of the optimal raw and class chain weights divided by the weighted key-set union of the query and target.
6. **Unknown-residue masking**: X, B, Z, J, O, `*`, and other nonstandard symbols are normalized to X. In the current main tracks, any k-mer covering X is excluded from anchor selection. Raw and Class share the same validity mask before class conversion.
7. **Reproducible indexes**: indexes record the actual parameters, positional postings, key frequencies, sequence IDF weights, algorithm version, and a BLAKE3 digest of the database.

## Installation and validation

Rust 1.85 or later is required, consistent with the minimum supported versions of the locked dependencies.

```bash
cargo install --path . --locked
cargo test
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

The release directory contains only the Rust implementation and requires neither Python nor Biopython.

## Quick start

```bash
# Synthetic-sequence demonstration; not biological validation evidence
chainmer demo

# Count unique seed keys for each sequence
chainmer seed -i db.faa -m chainmer

# Export TSV containing sequence IDs and 32-digit hexadecimal keys
chainmer seed -i db.faa -m subseqhash2 --keys -o subseqhash2.tsv

# Build an index and search using the parameters stored in that index
chainmer index -d db.faa -o db.chainmer.json --kc 5 --class-kc 8 --gmin 30 --gmax 300
chainmer search -d db.faa -q query.faa -x db.chainmer.json --top 5

# Compare seed counts and single-pass extraction times within the same Rust binary
chainmer compare -i db.faa -o seed_compare.tsv
```

Supported `seed --method` values are `chainmer`, `kmer`, `spaced`, `minimizer`, `syncmer`, `strobemer`, `ckmer`, and `subseqhash2`. See `src/cli.rs` for the fixed parameters of the conventional baselines. The syncmer baseline uses the **Closed** definition, in which the minimum s-mer occurs at the first or last position. A spaced pattern must contain at least one `1`. `simplified-strobemer` is a fixed-offset anchor pair, not a standard randstrobe or minstrobe.

By default, `search` outputs the new normalized collinear-chain score. `--raw` and `--class-only` run the respective single-track ablations. The optional weighted-search `--weights` JSON uses the legacy weighted-overlap score, with 32-digit hexadecimal keys as fields and finite floating-point numbers as values; it is not the new chaining model. `min_shared` must be at least 1, and zero-weight shared keys still produce neutral candidates. When an existing index is loaded, Chainmer validates the complete database contents, redundancy statistics, and cleaned-sequence coordinates rather than comparing only sequence IDs.

## Default k-mer parameters and tuning boundaries

The defaults `kc=5` and `class_kc=8` are **compromise operating points** for the current implementation, not a theoretically or empirically universal optimum. In the research repository, a 12-configuration full-factorial experiment evaluated raw `kc=4/5/6` and class `class_kc=6/7/8/9` on the same 177 [NiFe] hydrogenase sequences. Each configuration independently built an index and ran all queries three times. Raw 5 achieved the highest values of both candidate-precision metrics among the three raw lengths at every tested class length. No other configuration simultaneously dominated the default 5/8 setting across target-pair retention, both precision metrics, and candidate count. These results support retaining the current defaults, but they do not establish optimality across protein families.

With raw `kc=5` fixed, class length shows a clear retention–candidate-count trade-off:

| class_kc | Target-pair retention (%) | P_all (%) | P_Task3 (%) | Candidates/query |
|---|---:|---:|---:|---:|
| 6 | 88.48 | 15.84 | 71.97 | 69.073 |
| 7 | 80.35 | 22.67 | 74.11 | 43.819 |
| **8 (default)** | **74.31** | **33.48** | **80.98** | **27.435** |
| 9 | 67.55 | 39.46 | 83.69 | 21.164 |

Target-pair retention is based on 1,094 undirected ortholog pairs. P_all uses all unique candidate pairs as its denominator, including unlabeled pairs, and is therefore an operational candidate-precision metric. P_Task3 uses only candidate pairs with mutually exclusive ortholog/paralog labels. Candidate counts exclude self-matches. Consider class 7 when recall is prioritized, or class 9 when precision or candidate compression is more important. These are tuning directions observed on the current dataset only. New database deployments should be reevaluated using independent labels rather than selecting parameters from a single runtime measurement.

Changing `kc` or `class_kc` changes seed keys, DF/IDF, candidates, and index contents, so the index **must be rebuilt**. `search` uses the parameters stored in the index; changing only the search command is insufficient:

```bash
chainmer index -d db.faa -o db.kc5.class7.chainmer.json --kc 5 --class-kc 7
chainmer search -d db.faa -q query.faa -x db.kc5.class7.chainmer.json --top 5
```

The reproduction entry point in the research repository is `python OrthologySeed/scripts/rev_kc_sensitivity.py --repeats 3`. It writes `kc_sensitivity.tsv`, `kc_sensitivity_runs.tsv`, and `kc_sensitivity_metadata.json`. Benchmark sequences and labels are not distributed with the standalone Rust package. See Supplementary Data 5 of the paper for the complete experimental protocol and all 12 configurations.

## Rust API

Main entry points:

- `chainmer_seeds` / `chainmer_keys` / `ChainmerParams`: dual-track positional seed and key generation;
- `SeedIndex`: positional indexing, frequency statistics, saving, loading, and database validation;
- `search` / `CandidateScore`: IDF-normalized collinear-chain retrieval, per-track diagnostic scores, and positional-hit truncation statistics;
- `search_count` / `search_with_keys`: legacy inverted-index and exhaustive-count controls;
- `search_weighted_index` / `search_weighted_with_keys`: weighted retrieval;
- `contrastive_weights` / `specificity_weights` / `weighted_score`: training-fold weights and scoring;
- `subseqhash2_keys` / `SubseqHash2Params`: scalar protein adaptation of SubseqHash2.

All index queries check that query parameters match the index metadata. `gap_filter` and `gap_bucket` are independent switches; `set_use_gap` is only a shorthand for setting both.

## Scope of the SubseqHash2 adaptation

Based on the ABCk orders, forward and reverse dynamic programming, and pivot-merging concepts of the [official Shao-Group SubseqHash2 implementation](https://github.com/Shao-Group/subseqhash2), this package provides a **scalar Rust adaptation for the 20 standard amino acids plus X**. Its defaults are `window=25, subseq=11, d=5, repeats=11, table_seed=42, chunk_size=500`; these parameters were not tuned on the Chainmer test labels.

This is not the official four-character nucleotide SIMD binary, and its runtime cannot be directly equated with that of the official implementation. Windows containing X produce no seeds. Both forward and reverse DP have naive-subsequence enumeration controls, and chunked and whole-sequence processing have consistency tests.

## Hashing and compatibility

Hashes are stored as `u128`, but both 64-bit components derive from the same 64-bit rolling state. Therefore, the **effective state is at most 64 bits**; no 128-bit entropy or collision guarantee should be claimed. Anchor pairing uses only the lower 64 bits of each subhash and does not further reduce the information in the original state. The current index format is `chainmer-index/4`, and older indexes must be rebuilt. `chainmer_keys_legacy` retains only the v0.4 single-track key output without X masking for historical comparison; it is not used by the new index or default search.

FASTA input is parsed line by line with `BufRead`, avoiding an additional resident copy of the complete file string. The current API still returns a complete sequence map and is not a fully streaming proteome-scale index. Index saving first writes a temporary file in the same directory and performs `flush/sync_all`. Unix uses direct atomic replacement, while Windows uses backup and failure rollback to avoid losing an old index when replacement fails.

All 258,610 v0.4 Chainmer keys from the 177 benchmark sequences were checked individually, with zero sequences disagreeing between Rust and the Python prototype. This compatibility does not imply that v1.0.0 produces the same key sets or rankings as the older version.

## Pre-repair v0.5.0 feedback-subset smoke test (historical results)

The new default search was run on four `tree_input.fa` demonstration-tree subsets from the feedback package, excluding QUERY itself. All 7/7 archaeal references were retrieved with no additional candidates; 18/19 bacterial references were retrieved; all 17/17 fungal references appeared in the top 17; and the single viral reference was retrieved. The weak bacterial reference `WP_407854620.1` remained unretrieved. **These files are not the original full database. The results only validate execution and ranking behavior on real sequences; they are not full-database AP, FP/FN, or performance results.**

## Repaired benchmark summary (177-sequence reevaluation)

For the same 177 [NiFe] hydrogenase large-subunit sequences, genus-held-out evaluation, and deduplicated undirected pairs, repaired v1.0.0 computes database IDF from unlabeled whole-database statistics, while classification thresholds are still determined exclusively from training folds:

| Method | Task3 AP (5 folds) | Task4 AP (4 folds) |
|---|---:|---:|
| Chainmer v1.0.0 dual track | **0.825 ± 0.042** | **0.915 ± 0.056** |
| Chainmer v1.0.0 class only | 0.810 ± 0.034 | 0.903 ± 0.070 |
| Chainmer v1.0.0 raw only | 0.740 ± 0.041 | 0.863 ± 0.086 |

The dual-track AP values were 0.817 / 0.913 before repair and 0.825 / 0.915 in this reevaluation. This same-data comparison cannot isolate the contribution of individual repairs and does not replace independent cross-family validation. The dual track remains below the historical v0.4 SWPN results (0.916 / 0.949) and below BLASTp on Task4 (0.990). A single all-query run increased from 1.923 s to 4.546 s with strict chaining. See [BENCHMARK.md](BENCHMARK.md) for complete protocols, pre/post-repair tables, fold-level results, and limitations.

## License

MIT; see [LICENSE](LICENSE).
