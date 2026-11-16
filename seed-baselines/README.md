# seed-baselines — Rust implementations & benchmark data for all seed algorithms

Public repository: **https://github.com/rsj99-dev/Chainmer/tree/main/seed-baselines**

This directory contains the Rust implementations of **all eight seed algorithms** used in the paper
*"Using the Chainmer seed algorithm to help proteins find direct orthologs"* (Chainmer seed benchmark),
together with the **corresponding benchmark data**:

| Algorithm | Crate / module | CLI |
|---|---|---|
| k-mer (BLAST+) | root crate `src/seeds.rs::kmers` | `seed-baselines --method kmer --k 3 --stride 1` |
| k-mer (DIAMOND) | root crate `src/seeds.rs::spaced_kmers` | `seed-baselines --method spaced --masks 111101110111,111011010010111` |
| Minimizer | root crate `src/seeds.rs::minimizers` | `seed-baselines --method minimizer --k 7 --w 10` |
| Syncmer | root crate `src/seeds.rs::syncmers_with_mode` (closed) | `seed-baselines --method syncmer --k 7 --s 5` |
| Strobemer | root crate `src/seeds.rs::strobemers` (fixed-offset) | `seed-baselines --method strobemer --k 7 --w1 9 --w2 32 --t 2` |
| SubseqHash2 | root crate `src/subseqhash2.rs` (protein scalar) | `seed-baselines --method subseqhash2 --window 25 --subseq 11 --d 5 --repeats 11` |
| **Chainmer** | `chainmer/` (full package v1.0.2, δ=5) | `chainmer seed --method chainmer --keys …` |
| **Chainmer-op1** | `chainmer-op1/` (full package v1.0.0, five-group class channel, δ=3) | `chainmer-op1 seed --method chainmer --keys …` |

## Layout

```
seed-baselines/
├── Cargo.toml, src/, tests/      six-baseline crate (lib seed_baselines + bin seed-baselines),
│                                 byte-exact compatible with the paper exporter; 6 golden regression tests
├── chainmer/                     Chainmer v1.0.2 full package (dual-channel anchor-pair seeds, index/search)
├── chainmer-op1/                 Chainmer-op1 v1.0.0 full package (five-group full-coverage class channel, δ=3)
├── data/
│   ├── hydrogenase_datasets/     benchmark INPUT data: HydDB 2012, genus-level hold-out, KEGG 92
│   │                             (136 files, SHA-256 MANIFEST.tsv)
│   └── results/                  benchmark OUTPUT data for all 8 methods: metrics_toolB/native, analysis,
│                                 summary CSVs, full2012 FP/FN, McNemar, fig5 tree cases, per-query preds,
│                                 gene_tree_rooted.nwk, reproduction_check.md
├── README.md / README_CN.md / CHANGELOG.md / LICENSE / .gitignore
```

## Why Rust, and why no SIMD

- Cross-platform compatibility: one Rust codebase, consistent semantics/compilation across OS/arch; no C++
  undefined behavior or compiler/template divergence.
- Deterministic computation: x86 CPU vector acceleration deliberately disabled — pure scalar portable code,
  no SIMD intrinsics (SubseqHash2 is a protein scalar adaptation). Key sets and scores are byte-reproducible.
- Engineering reliability: compile-time borrow/lifetime checks eliminate memory-class bugs before runtime.
- Experimental consistency: Chainmer itself is Rust; all baselines share the same hash (splitmix64), cleaning
  (U→C, non-standard→X), dedup and scoring code, minimizing implementation-layer error.

## Data notes

- `data/hydrogenase_datasets/` is the exact input data used in the paper (SHA-256 verified).
- `data/results/` contains the paper's benchmark outputs for **all eight methods, including Chainmer and
  Chainmer-op1** (accuracy, Macro-F1, MRR, Recall@5, AUC, answer rate, confusion matrices, per-group recall,
  full-database FP/FN, McNemar tests, fig5 tree cases, per-query predictions).
- Per-sequence seed-key TSVs are GB-scale and are **not** stored here; they are regenerable **byte-exactly**
  via the CLIs above (verified 57/57 SHA-256 identical to the paper's in-use exports; see
  `data/results/reproduction_check.md`).

## Build & test

```
cargo build --release            # six-baseline crate (root)
cargo test  --release            # 6 golden regression tests + unit tests
cd chainmer      && cargo build --release   # Chainmer v1.0.2
cd chainmer-op1  && cargo build --release   # Chainmer-op1 v1.0.0
```

On Chinese-character paths use the rust-lld workaround (link-self-contained) as documented in each crate.

## License

MIT (see LICENSE).
