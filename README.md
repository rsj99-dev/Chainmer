<p align="center">
  <img src="Chainmer_Logo.png" alt="Chainmer Logo" width="420">
</p>

# Chainmer v0.5.1

Chainmer is a Rust library and command-line tool for **protein homolog candidate retrieval and prefiltering**. v0.5.1 replaces the previous “single-channel class-equivalent key counting” with **raw/class dual-channel positional seeds, database IDF calibration, and collinear chain scoring**; the old counting and SWPN interfaces are retained for ablation comparisons.

> A Seed Hit only indicates that candidates share seeds; it does not constitute a homolog or ortholog determination. Orthology requires downstream analyses such as gene tree–species tree reconciliation; HydDB subgroup labels are also not equivalent to ortholog groups.

## Core Components

1. **Dual-channel recall**: the raw channel uses original amino acids with `kc=5`, and the class channel uses substitution equivalence classes with `class_kc=8`; the two kinds of keys use different hash salts;
2. **Positional anchor pairs**: stores dual anchor positions for the query and target, and no longer treats seeds merely as unordered key sets;
3. **Database frequency calibration**: computes IDF from the current database document frequency, and automatically downweights high-frequency keys;
4. **Duplicate-contribution constraint**: retained positional hits for the same key within a candidate jointly share one IDF weight; each key retains at most 64 representative hits via stable-hash bottom-k, the retained set does not depend on posting order, and counts before and after truncation are reported;
5. **Collinear chain scoring**: within staggered diagonal buckets, finds the optimal chain under the compatibility condition that all four anchor endpoints are strictly increasing in both query and target; the final score is “sum of raw/class optimal chain weights ÷ union of weighted key sets of query and target”;
6. **Unknown residue masking**: X, B, Z, J, O, `*`, and other non-standard symbols are sanitized to X; the current main channels do not allow any k-mer covering X to participate in anchor selection, and Raw/Class share the same validity mask before class conversion;
7. **Reproducible index**: records actual parameters, positional postings, key frequencies, sequence IDF weights, algorithm version, and database BLAKE3 digest.

## Installation and Verification

Requires Rust 1.85 or higher (consistent with the minimum requirements of the locked dependencies).

```bash
cargo install --path . --locked
cargo test
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

The release directory contains only the Rust implementation and does not require Python or Biopython.

## Quick Start

```bash
# Synthetic sequence demo; not biological validation evidence
chainmer demo

# Count unique seed keys per sequence
chainmer seed -i db.faa -m chainmer

# Export TSV: sequence ID + 32-bit hexadecimal key
chainmer seed -i db.faa -m subseqhash2 --keys -o subseqhash2.tsv

# Build an index and search using parameters saved in the index
chainmer index -d db.faa -o db.chainmer.json --kc 5 --class-kc 8 --gmin 30 --gmax 300
chainmer search -d db.faa -q query.faa -x db.chainmer.json --top 5

# Compare seed counts and single-extraction time within the same Rust binary
chainmer compare -i db.faa -o seed_compare.tsv
```

Supported `seed --method`: `chainmer`, `kmer`, `spaced`, `minimizer`, `syncmer`, `strobemer`, `ckmer`, `subseqhash2`. Fixed parameters for the traditional baselines are in `src/cli.rs`; the syncmer baseline uses the **Closed** definition where the minimum s-mer is at the first or last end, and a spaced pattern must contain at least one `1`; `simplified-strobemer` is a fixed-offset anchor-pair, not a standard randstrobe/minstrobe.

Default `search` outputs the new normalized collinear-chain score. `--raw` and `--class-only` perform single-channel ablations, respectively. The `--weights` JSON for weighted retrieval follows the legacy weighted-overlap scoring, using 32-bit hexadecimal keys as field names and finite floats as values; it is not the new chaining model. `min_shared` must be at least 1; zero-weight shared keys still form neutral candidates. When loading an existing index, the database’s full contents, redundant statistics, and sanitized sequence coordinates are validated, rather than only comparing sequence IDs.

## Rust API

Main entry points:

- `chainmer_seeds` / `chainmer_keys` / `ChainmerParams`: dual-channel positional seed and key generation;
- `SeedIndex`: positional index, frequency statistics, saving, loading, and database validation;
- `search` / `CandidateScore`: IDF-normalized collinear chain retrieval, per-channel diagnostic scores, and positional hit truncation statistics;
- `search_count` / `search_with_keys`: legacy inverted-index and exhaustive counting comparisons;
- `search_weighted_index` / `search_weighted_with_keys`: weighted retrieval;
- `contrastive_weights` / `specificity_weights` / `weighted_score`: training-fold weights and scoring;
- `subseqhash2_keys` / `SubseqHash2Params`: SubseqHash2 protein scalar adaptation.

All index queries check whether query parameters are consistent with index metadata. `gap_filter` and `gap_bucket` are independent switches; `set_use_gap` is only a shorthand that sets both at once.

## SubseqHash2 Adaptation Boundary

This package provides a **scalar Rust adaptation for the 20 standard amino acids + X**, based on the ABCk orders, forward/reverse dynamic programming, and pivot merging ideas of the [Shao-Group official SubseqHash2 implementation](https://github.com/Shao-Group/subseqhash2). The default parameters are `window=25, subseq=11, d=5, repeats=11, table_seed=42, chunk_size=500`, and have not been tuned on Chainmer test labels.

It is not the official four-character nucleotide SIMD binary and cannot be directly equated with the official implementation’s runtime. Windows containing X do not produce seeds. Both forward and reverse DP have naive subsequence enumeration comparison tests, and chunked vs. whole-segment processing has consistency tests.

## Hash and Compatibility

The hash is stored as `u128`, but the two 64-bit components come from the same 64-bit rolling state, so **the effective state is at most 64 bit**; it must not be claimed to have 128-bit entropy or collision guarantees; anchor combinations only use the lower 64 bits of the sub-hash and do not further reduce the information content of the original state. The current index format is `chainmer-index/4`, and old indexes must be rebuilt. `chainmer_keys_legacy` retains only v0.4 single-channel key output with X unmasked for historical checking and does not participate in new indexes or default retrieval.

FASTA uses `BufRead` line-by-line parsing to avoid additionally keeping the full file string resident; the current API still returns the full sequence map and is not a fully streaming proteome index. Index saving first writes a temporary file in the same directory and performs `flush/sync_all`; Unix uses direct atomic replacement, while Windows uses backup and rollback-on-failure to avoid losing the old index if overwriting fails.

177 benchmark sequences and 258,610 v0.4 Chainmer keys have been checked item by item; the number of sequences where the Rust and Python prototypes disagree is 0; this compatibility does not mean that v0.5.1 produces the same key sets or rankings as the old version.

## Pre-fix v0.5.0 Feedback-Subset Smoke Test (Historical Results)

Using the four `tree_input.fa` display-tree subsets from the feedback package, the new default retrieval was actually run, excluding QUERY itself: archaea 7/7 references all entered the candidates and no other candidates were reported; bacteria 18/19; fungi 17/17 within the top 17; viruses 1/1. The weak bacterial reference `WP_407854620.1` was still not recalled. **These files are not the original full database; the results only verify runtime and ranking behavior on real sequences, not full-database AP, FP/FN, or performance conclusions.**

## Fixed-Version Benchmark Summary (177-Sequence Re-evaluation)

Under the same 177 [NiFe] hydrogenase large-subunit sequences, leave-genus-out labels, and undirected deduplicated pairs, the database IDF of fixed v0.5.1 uses full-database unlabeled statistics, and the classification threshold is still determined only by the training folds:

| Method | Task3 AP (5-fold) | Task4 AP (4-fold) |
|---|---:|---:|
| Chainmer v0.5.1 dual-channel | **0.825 ± 0.042** | **0.915 ± 0.056** |
| Chainmer v0.5.1 class single-channel | 0.810 ± 0.034 | 0.903 ± 0.070 |
| Chainmer v0.5.1 raw single-channel | 0.740 ± 0.041 | 0.863 ± 0.086 |

Before the fix, the dual-channel results were 0.817 / 0.913; this round’s re-evaluation gives 0.825 / 0.915; this same-data comparison cannot decompose the contribution of individual fixes, nor does it replace independent cross-family validation. The dual-channel result still does not exceed the historical v0.4 SWPN result (0.916 / 0.949) or BLASTp on Task4 (0.990). After stricter chaining, the single-run time for all queries increased from 1.923 s to 4.546 s. For full definitions, pre-/post-fix tables, per-fold results, and limitations, see [BENCHMARK.md](BENCHMARK.md).

## License

MIT, see [LICENSE](LICENSE).
```
