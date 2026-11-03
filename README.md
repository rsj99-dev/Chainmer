# EvoSeed v0.5.0

EvoSeed 是一个用于**蛋白同源/直系同源候选召回与排序**的 Rust 库和命令行工具。v0.5 将原来的“单通道类等价键计数”替换为 **raw/class 双通道位置种子、数据库 IDF 校准和共线链评分**；旧计数与 SWPN 接口保留作消融对照。

> Seed Hit 只表示候选共享种子，不构成同源或直系同源判定。直系同源需要基因树—物种树协调等后续分析；HydDB 亚群标签也不等同于直系同源群。

## 核心组成

1. **双通道召回**：raw 通道使用原始氨基酸 `kc=5`，class 通道使用替换等价类 `class_kc=8`，两种键使用不同哈希盐；
2. **位置化锚点对**：保存查询与目标的双锚点位置，不再把种子仅视为无序键集合；
3. **数据库频率校准**：以当前数据库 document frequency 计算 IDF，高频键自动降权；
4. **重复贡献约束**：同一键在一个候选中的全部位置命中分摊一份 IDF 权重，并限制每键位置配对数；
5. **共线链评分**：在错开的对角线分桶内求单调最优链，最终分数为“raw/class 最优链权重之和 ÷ 查询与目标的加权键集并集”；
6. **可复现索引**：记录实际参数、位置 postings、键频率、序列 IDF 权重、算法版本和数据库 BLAKE3 摘要。

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
evoseed index -d db.faa -o db.evoseed.json --kc 5 --class-kc 8 --gmin 30 --gmax 300
evoseed search -d db.faa -q query.faa -x db.evoseed.json --top 5

# 同一 Rust 二进制内的种子数量与单次提取耗时比较
evoseed compare -i db.faa -o seed_compare.tsv
```

支持的 `seed --method`：`evoseed`、`kmer`、`spaced`、`minimizer`、`syncmer`、`strobemer`、`ckmer`、`subseqhash2`。传统基线的固定参数见 `src/cli.rs`；`simplified-strobemer` 是 fixed-offset anchor-pair，不是标准 randstrobe/minstrobe。

默认 `search` 输出新的共线链归一化分数。`--raw` 和 `--class-only` 分别执行单通道消融。加权检索的 `--weights` JSON 沿用旧式加权重叠评分，以 32 位十六进制键为字段名、有限浮点数为值；它不是新的链化模型。`min_shared` 必须至少为 1；零权重共享键仍形成中性候选。加载已有索引时会校验数据库完整内容，而不是只比较序列 ID。

## Rust API

主要入口：

- `evoseed_seeds` / `evoseed_keys` / `EvoParams`：双通道位置种子与键生成；
- `SeedIndex`：位置索引、频率统计、保存、加载和数据库校验；
- `search` / `CandidateScore`：IDF 归一化的共线链检索及分通道诊断分；
- `search_count` / `search_with_keys`：旧式倒排与穷举计数对照；
- `search_weighted_index` / `search_weighted_with_keys`：加权检索；
- `contrastive_weights` / `specificity_weights` / `weighted_score`：训练折权重与评分；
- `subseqhash2_keys` / `SubseqHash2Params`：SubseqHash2 蛋白标量适配。

所有索引查询都会检查查询参数与索引元数据是否一致。`gap_filter` 与 `gap_bucket` 是独立开关；`set_use_gap` 只是同时设置两者的简写。

## SubseqHash2 适配边界

本包依据 [Shao-Group 官方 SubseqHash2 实现](https://github.com/Shao-Group/subseqhash2)的 ABCk orders、前向/反向动态规划和 pivot 合并思路，提供 **20 种标准氨基酸 + X 的标量 Rust 适配**。默认参数为 `window=25, subseq=11, d=5, repeats=11, table_seed=42, chunk_size=500`，未在 EvoSeed 测试标签上调优。

它不是官方四字符核酸 SIMD 二进制，不能与官方实现的运行时间直接等同。包含 X 的窗口不产种子。前向和反向 DP 均有朴素子序列枚举对照测试，分块与整段处理有一致性测试。

## 哈希与兼容性

哈希以 `u128` 存储，但两个 64 位分量来自同一个 64 位滚动状态，因此**有效状态最多 64 bit**，不得宣称具有 128-bit 熵或碰撞保证。v0.5 索引格式为 `evoseed-index/3`，旧版索引必须重建。`evoseed_keys_legacy` 仅保留 v0.4 单通道键输出用于历史核对，不参与新索引和默认检索。

177 条基准序列、258,610 个 v0.4 EvoSeed 键已逐项核对，Rust 与 Python 原型不一致的序列数为 0；该兼容性不代表 v0.5 与旧版产生相同键集或排名。

## v0.5 反馈子集 smoke test

使用反馈包四个 `tree_input.fa` 展示树子集实际运行新默认检索，排除 QUERY 自身后：古菌 7/7 参照全部进入候选且未报出其他候选；细菌 18/19；真菌 17/17 位于前 17；病毒 1/1。细菌弱参照 `WP_407854620.1` 仍未召回。**这些文件不是原始全库，结果只验证真实序列上的运行和排序行为，不是全库 AP、FP/FN 或性能结论。**

## 基准摘要

以下数值均为 **v0.4 历史结果**，尚未用 v0.5 双通道链化算法重跑，不得用于宣称新算法性能。

同一 177 条 [NiFe] 氢化酶大亚基数据、留属划分、无向去重配对和训练折阈值流程下：

| 方法 | Task3 AUPRC（5 折） | Task4 AUPRC（4 折） |
|---|---:|---:|
| EvoSeed-SWPN | 0.916 ± 0.045 | 0.949 ± 0.045 |
| EvoSeed 计数 | 0.802 ± 0.059 | 0.919 ± 0.050 |
| SubseqHash2 蛋白标量适配 | 0.719 ± 0.043 | 0.867 ± 0.087 |

AUPRC 使用按同分组的非插值平均精度（AP），AUROC 同样按同分组计算。此数据上 v0.4 EvoSeed-SWPN 在 Task3 最高，BLASTp 在 Task4 最高（0.990）；SubseqHash2 适配在两个任务上均未超过 v0.4 EvoSeed 计数或 EvoSeed-SWPN。该结果不支持脱离数据集的普遍排名，也不代表 v0.5。

本次 Rust release 单次顺序提取中，EvoSeed 为 1461.1 键/序列、29.3 ms/177 序列；SubseqHash2 为 2776.3 键/序列、5388.0 ms。后者候选数较少，但键数和标量提取成本更高。完整口径、命令与限制见 [BENCHMARK.md](BENCHMARK.md)。

## 许可

MIT，见 [LICENSE](LICENSE)。
