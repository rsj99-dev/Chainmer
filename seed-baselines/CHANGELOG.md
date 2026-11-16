# Changelog

## v1.0.0 (2026-11-16)

- 首发：六类基线种子算法（kmer / spaced / minimizer / syncmer / strobemer / subseqhash2）
  标量 Rust 实现，与论文 tools/seed_dump.rs 逐字节兼容。
- 不含 x86 SIMD intrinsics（纯标量可移植），跨平台键集合确定性可复现。
- 六项回归 golden 测试（论文在用二进制生成）＋清洗/参数校验单元测试。
- 共享 Chainmer 内核的 splitmix64 哈希、U→C/非标准→X 清洗与 BTreeSet 去重语义。
