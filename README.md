# EvoSeed v0.4.0

EvoSeed 是一个用于**蛋白同源/直系同源候选召回与排序**的 Rust 库和命令行工具。它组合替换等价类、类 minimizer 锚点、有序锚点对、间距过滤与分桶，并支持仅由训练数据拟合的 SW/SWPN 权重。

> Seed Hit 只表示候选共享种子，不构成同源或直系同源判定。直系同源需要基因树—物种树协调等后续分析；HydDB 亚群标签也不等同于直系同源群。

## 核心组成

1. **替换等价类**：20 种标准氨基酸分入 5 个互斥类，特殊或未知残基进入第 6 类；
2. **有序锚点对**：在类串或原始氨基酸串上选择 minimizer 锚点，组合有序锚点对；
3. **间距约束**：可分别控制距离过滤和间距分桶；
4. **对比加权**：SWPN 为 `ln((f_pos + 1) / (f_neg + 1))`，未出现的键默认权重为 0；
5. **可复现索引**：记录实际参数、类别映射、算法版本、哈希方案和数据库内容 BLAKE3 摘要。

## 安装与验证

需要 Rust 1.85 或更高版本（与锁定依赖的最低要求一致）。

```bash
cargo install --path . --locked
cargo test
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

发布目录只包含 Rust 实现，不需要 Python 或 Biopython。

## 快速开始

```bash
# 合成序列演示，不是生物学验证证据
evoseed demo

# 统计每条序列的唯一种子键
evoseed seed -i db.faa -m evoseed

# 导出 TSV：序列 ID + 32 位十六进制键
evoseed seed -i db.faa -m subseqhash2 --keys -o subseqhash2.tsv

# 构建索引并使用索引中保存的参数检索
evoseed index -d db.faa -o db.evoseed.json --kc 5 --gmin 30 --gmax 300
evoseed search -d db.faa -q query.faa -x db.evoseed.json --top 5

# 同一 Rust 二进制内的种子数量与单次提取耗时比较
evoseed compare -i db.faa -o seed_compare.tsv
```

支持的 `seed --method`：`evoseed`、`kmer`、`spaced`、`minimizer`、`syncmer`、`strobemer`、`ckmer`、`subseqhash2`。传统基线的固定参数见 `src/cli.rs`；`simplified-strobemer` 是 fixed-offset anchor-pair，不是标准 randstrobe/minstrobe。

加权检索的 `--weights` JSON 以 32 位十六进制键为字段名、有限浮点数为值。`min_shared` 必须至少为 1；零权重共享键仍形成中性候选。加载已有索引时会校验数据库完整内容，而不是只比较序列 ID。

## Rust API

主要入口：

- `evoseed_keys` / `anchors` / `EvoParams`：EvoSeed 键生成；
- `SeedIndex`：索引构建、保存、加载和数据库校验；
- `search` / `search_with_keys`：倒排与穷举计数检索；
- `search_weighted_index` / `search_weighted_with_keys`：加权检索；
- `contrastive_weights` / `specificity_weights` / `weighted_score`：训练折权重与评分；
- `subseqhash2_keys` / `SubseqHash2Params`：SubseqHash2 蛋白标量适配。

所有索引查询都会检查查询参数与索引元数据是否一致。`gap_filter` 与 `gap_bucket` 是独立开关；`set_use_gap` 只是同时设置两者的简写。

## SubseqHash2 适配边界

本包依据 [Shao-Group 官方 SubseqHash2 实现](https://github.com/Shao-Group/subseqhash2)的 ABCk orders、前向/反向动态规划和 pivot 合并思路，提供 **20 种标准氨基酸 + X 的标量 Rust 适配**。默认参数为 `window=25, subseq=11, d=5, repeats=11, table_seed=42, chunk_size=500`，未在 EvoSeed 测试标签上调优。

它不是官方四字符核酸 SIMD 二进制，不能与官方实现的运行时间直接等同。包含 X 的窗口不产种子。前向和反向 DP 均有朴素子序列枚举对照测试，分块与整段处理有一致性测试。

## 哈希与兼容性

为保持已有 EvoSeed 基准可比，v0.4.0 保留原型的键输出。哈希以 `u128` 存储，但两个 64 位分量来自同一个 64 位滚动状态，因此**有效状态最多 64 bit**，不得宣称具有 128-bit 熵或碰撞保证。索引格式为 `evoseed-index/2`，不兼容旧版无完整元数据的索引。

177 条基准序列、258,610 个 EvoSeed 键已逐项核对，Rust 与 Python 原型不一致的序列数为 0。

## 基准摘要

同一 177 条 [NiFe] 氢化酶大亚基数据、留属划分、无向去重配对和训练折阈值流程下：

| 方法 | Task3 AUPRC（5 折） | Task4 AUPRC（4 折） |
|---|---:|---:|
| EvoSeed-SWPN | 0.916 ± 0.045 | 0.949 ± 0.045 |
| EvoSeed 计数 | 0.802 ± 0.059 | 0.919 ± 0.050 |
| SubseqHash2 蛋白标量适配 | 0.719 ± 0.043 | 0.867 ± 0.087 |

AUPRC 使用按同分组的非插值平均精度（AP），AUROC 同样按同分组计算。此数据上 EvoSeed-SWPN 在 Task3 最高，BLASTp 在 Task4 最高（0.990）；SubseqHash2 适配在两个任务上均未超过 EvoSeed 计数或 EvoSeed-SWPN。该结果不支持脱离数据集的普遍排名。

本次 Rust release 单次顺序提取中，EvoSeed 为 1461.1 键/序列、29.3 ms/177 序列；SubseqHash2 为 2776.3 键/序列、5388.0 ms。后者候选数较少，但键数和标量提取成本更高。完整口径、命令与限制见 [BENCHMARK.md](BENCHMARK.md)。

## 许可

MIT，见 [LICENSE](LICENSE)。
