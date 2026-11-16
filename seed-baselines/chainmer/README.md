<p align="center">
  <img src="Chainmer_Logo.png" alt="Chainmer Logo" width="420">
</p>

# Chainmer v1.0.2

**Language:** English | [中文](README_CN.md)

Chainmer is a Rust library and command-line tool for **protein homology candidate retrieval and prefiltering**. v1.0.0 replaced the former single-track class-equivalence key count with **raw/class dual-track positional seeds, database-calibrated IDF, and collinear chain scoring**. **v1.0.2 changes the default gap-bin width `delta` from 15 to 5** (3× finer gap binning, fewer chance paralog matches); the rest of the pipeline is unchanged from v1.0.1. Legacy count and SWPN interfaces remain as ablation controls.

> A Seed Hit only indicates that candidates share a seed; it does not establish homology or orthology. Orthology requires downstream analyses such as gene tree–species tree reconciliation. HydDB subgroup labels are not equivalent to orthologous groups.

## Core components

1. **Dual-track recall**: raw track uses original amino acids with `kc=5`; class track uses substitution-equivalence classes with `class_kc=8`; the two key types use different hash salts.
2. **Positioned anchor pairs**: both anchor positions are stored for queries and targets.
3. **Database frequency calibration**: IDF from current-database document frequencies; frequent keys down-weighted.
4. **Duplicate contribution constraint**: retained positional hits of one key within a candidate share a single IDF weight; stable-hash bottom-k keeps ≤64 representative hits per key.
5. **Collinear chain scoring**: optimal chains within staggered diagonal bins, requiring all four anchor endpoints to increase strictly in both query and target.
6. **Unknown-residue masking**: X, B, Z, J, O, `*`, and other nonstandard symbols normalized to X; k-mers covering X are excluded from anchor selection.
7. **Reproducible indexes**: record actual parameters, positional postings, key frequencies, sequence IDF weights, algorithm version, and a BLAKE3 database digest.

## Installation and validation

Rust 1.85 or later required.

```bash
cargo install --path . --locked
cargo test
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

## Quick start

```bash
chainmer demo                                   # synthetic demo, not biological validation
chainmer seed -i db.faa -m chainmer             # per-sequence unique seed-key counts
chainmer seed -i db.faa -m chainmer --keys -o chainmer.tsv   # TSV export (paper key protocol)
chainmer index -d db.faa -o db.chainmer.idx --kc 5 --class-kc 8 --gmin 30 --gmax 300 --delta 5
chainmer search -d db.faa -q query.faa -x db.chainmer.idx --top 5
chainmer compare -i db.faa -o seed_compare.tsv  # seed counts / extraction time within one binary
```

Supported `seed --method`: `chainmer`, `kmer`, `spaced`, `minimizer`, `syncmer`, `strobemer`, `ckmer`, `subseqhash2`. The syncmer baseline uses the **Closed** definition; spaced patterns must contain at least one `1`; `simplified-strobemer` is a fixed-offset anchor pair (not standard randstrobe/minstrobe).

## v1.0.2 defaults

`kc=5`, `class_kc=8`, `w1=6`, gap window `[30,300]`, **`delta=5`**, ≤12 pairs per anchor. `delta` changed from 15 (v1.0.1) to 5 — the only configuration weakly Pareto-dominating the old default in a 30-configuration three-stage scan (see CHANGELOG v1.0.2). Changing `kc`/`class_kc`/`delta` changes seed keys and index contents: **rebuild the index**.

## v1.0.2 benchmark (HydDB 2012 / genus hold-out / KEGG 92, seed-key retrieval protocol)

Protocol: unified IDF-weighted Jaccard over seed-key sets (the paper's main protocol), leave-one-out / genus-level hold-out. Data and all-method outputs: [seed-baselines/data](https://github.com/rsj99-dev/Chainmer/tree/main/seed-baselines) (`data/hydrogenase_datasets/`, `data/results/`).

| Task | Top-1 (95% CI) | Macro-F1 | MRR | Recall@5 | AUC | Answer rate |
|---|---:|---:|---:|---:|---:|---:|
| HydDB 2012 full | **0.9945** (0.991–0.998) | 0.9896 | 0.9954 | 0.9955 | 0.9244 | 1.000 |
| Genus hold-out | **0.9275** (0.900–0.950) | 0.9013 | 0.9476 | 0.9675 | 0.8896 | 1.000 |
| KEGG 92 external | **1.0000** | 1.0000 | 1.0000 | 1.0000 | 0.8478 | 1.000 |

Full-2012 per-query paralog analysis: **false positives (paralog over-detection) 295 (14.7%), false negatives (ortholog miss) 0**; median density **6.049 keys/residue**, median 3086 keys/sequence.

Hold-out McNemar (Chainmer vs others): vs k-mer (BLAST+) 16:18 (p=0.864, n.s.); vs Minimizer 43:6, Syncmer 41:17 (p=0.0022), SubseqHash2 53:11, Strobemer 174:2, k-mer (DIAMOND) 40:14 (p=0.0005) all significant. Chainmer significantly exceeds all general-purpose seeds except exhaustive k-mer (BLAST+).

> Note: the table above is the **seed-key retrieval** protocol (paper main comparison). This package's `search` collinear-chain score is a different protocol; the v1.0.0/v1.0.1 177-sequence collinear-chain AP historical results are in [BENCHMARK.md](BENCHMARK.md) (historical appendix) and CHANGELOG and must not be conflated with the table above.

## Rust API

- `chainmer_seeds` / `chainmer_keys` / `ChainmerParams`: dual-track positional seeds and key generation;
- `SeedIndex`: positional indexing, frequency statistics, save/load, database validation;
- `search` / `CandidateScore`: IDF-normalized collinear-chain retrieval, per-track diagnostics, hit-truncation statistics;
- `search_count` / `search_with_keys`: legacy inverted-index and exhaustive-count controls;
- `search_weighted_index` / `search_weighted_with_keys`: weighted retrieval;
- `subseqhash2_keys` / `SubseqHash2Params`: scalar protein adaptation of SubseqHash2.

## Scope of the SubseqHash2 adaptation

Based on the ABCk orders, forward/reverse DP, and pivot-merging concepts of the [official Shao-Group implementation](https://github.com/Shao-Group/subseqhash2), this package provides a **scalar Rust adaptation for the 20 standard amino acids plus X** (defaults `window=25, subseq=11, d=5, repeats=11, table_seed=42, chunk_size=500`). It is not the official four-character nucleotide SIMD binary; windows containing X produce no seeds.

## Hashing and compatibility

Hashes are stored as `u128` but the effective state is at most 64 bits; no 128-bit entropy claim. Current index format `chainmer-index/5` (compact binary CSR); older indexes must be rebuilt. `chainmer_keys_legacy` retains only the v0.4 single-track key output for historical comparison.

## License

MIT; see [LICENSE](LICENSE).
