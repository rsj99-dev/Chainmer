//! seed-baselines：论文六类基线种子算法的标量、无 SIMD 可移植 Rust 实现。
//!
//! 与论文导出工具（paper_seed_benchmark/tools/seed_dump.rs）逐字节兼容：
//! 同一哈希（splitmix64）、清洗（U→C、非标准→X）、去重与参数语义。
//! 不含 x86 向量加速指令集（纯标量），保证跨平台键集合逐字节可复现。

pub mod fasta;
pub mod seeds;
pub mod subseqhash2;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
