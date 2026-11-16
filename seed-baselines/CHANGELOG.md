# Changelog

## repo v1.1.0 (2026-11-16)

- 仓库扩展为全部八类种子算法的单仓库：新增 chainmer/（v1.0.2 完整包）与 chainmer-op1/（v1.0.0 完整包）子 crate。
- 新增 data/hydrogenase_datasets/（基准输入，136 文件 SHA-256）与 data/results/（全部 8 方法基准输出，含 Chainmer/op1）。
- 公开地址：https://github.com/rsj99-dev/Chainmer/tree/main/seed-baselines 。
- 六基线 crate 保持 v1.0.0，逐字节兼容与 6 项黄金标准回归不变。

## v1.0.0 (2026-11-16)

- 首发：六类基线种子算法（kmer / spaced / minimizer / syncmer / strobemer / subseqhash2）
  标量 Rust 实现，与论文 tools/seed_dump.rs 逐字节兼容。
- 不含 x86 SIMD intrinsics（纯标量可移植），跨平台键集合确定性可复现。
- 六项回归 golden 测试（论文在用二进制生成）＋清洗/参数校验单元测试。
- 共享 Chainmer 内核的 splitmix64 哈希、U→C/非标准→X 清洗与 BTreeSet 去重语义。
