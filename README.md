<p align="center">
  <img src="Chainmer_Logo.png" alt="Chainmer Logo" width="420">
</p>

# Chainmer v0.5.1

Chainmer 是一个用于**蛋白同源候选检索与预过滤**的 Rust 库和命令行工具。v0.5.1 将原来的“单通道类等价键计数”替换为 **raw/class 双通道位置种子、数据库 IDF 校准和共线链评分**；旧计数与 SWPN 接口保留作消融对照。

> Seed Hit 只表示候选共享种子，不构成同源或直系同源判定。直系同源需要基因树—物种树协调等后续分析；HydDB 亚群标签也不等同于直系同源群。

## 核心组成

1. **双通道召回**：raw 通道使用原始氨基酸 `kc=5`，class 通道使用替换等价类 `class_kc=8`，两种键使用不同哈希盐；
2. **位置化锚点对**：保存查询与目标的双锚点位置，不再把种子仅视为无序键集合；
3. **数据库频率校准**：以当前数据库 document frequency 计算 IDF，高频键自动降权；
4. **重复贡献约束**：同一键在一个候选中的保留位置命中共同分摊一份 IDF 权重；每键以稳定哈希 bottom-k 最多保留 64 个代表命中，保留集合不依赖 posting 顺序，并报告截断前后数量；
5. **共线链评分**：在错开的对角线分桶内，以四个锚点端点在查询和目标中均严格递增为兼容条件求最优链，最终分数为“raw/class 最优链权重之和 ÷ 查询与目标的加权键集并集”；
6. **未知残基屏蔽**：X、B、Z、J、O、`*` 和其他非标准符号清洗为 X；当前主通道不允许任何覆盖 X 的 k-mer 参与锚点选择，Raw/Class 在类别转换前共用同一有效性掩码；
7. **可复现索引**：记录实际参数、位置 postings、键频率、序列 IDF 权重、算法版本和数据库 BLAKE3 摘要。

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
chainmer demo

# 统计每条序列的唯一种子键
chainmer seed -i db.faa -m chainmer

# 导出 TSV：序列 ID + 32 位十六进制键
chainmer seed -i db.faa -m subseqhash2 --keys -o subseqhash2.tsv

# 构建索引并使用索引中保存的参数检索
chainmer index -d db.faa -o db.chainmer.json --kc 5 --class-kc 8 --gmin 30 --gmax 300
chainmer search -d db.faa -q query.faa -x db.chainmer.json --top 5

# 同一 Rust 二进制内的种子数量与单次提取耗时比较
chainmer compare -i db.faa -o seed_compare.tsv
```

支持的 `seed --method`：`chainmer`、`kmer`、`spaced`、`minimizer`、`syncmer`、`strobemer`、`ckmer`、`subseqhash2`。传统基线的固定参数见 `src/cli.rs`；syncmer 基线使用最小 s-mer 位于首端或末端的 **Closed** 定义，spaced pattern 至少要包含一个 `1`；`simplified-strobemer` 是 fixed-offset anchor-pair，不是标准 randstrobe/minstrobe。

默认 `search` 输出新的共线链归一化分数。`--raw` 和 `--class-only` 分别执行单通道消融。加权检索的 `--weights` JSON 沿用旧式加权重叠评分，以 32 位十六进制键为字段名、有限浮点数为值；它不是新的链化模型。`min_shared` 必须至少为 1；零权重共享键仍形成中性候选。加载已有索引时会校验数据库完整内容、冗余统计及清洗后序列坐标，而不是只比较序列 ID。

## Rust API

主要入口：

- `chainmer_seeds` / `chainmer_keys` / `ChainmerParams`：双通道位置种子与键生成；
- `SeedIndex`：位置索引、频率统计、保存、加载和数据库校验；
- `search` / `CandidateScore`：IDF 归一化的共线链检索、分通道诊断分及位置命中截断统计；
- `search_count` / `search_with_keys`：旧式倒排与穷举计数对照；
- `search_weighted_index` / `search_weighted_with_keys`：加权检索；
- `contrastive_weights` / `specificity_weights` / `weighted_score`：训练折权重与评分；
- `subseqhash2_keys` / `SubseqHash2Params`：SubseqHash2 蛋白标量适配。

所有索引查询都会检查查询参数与索引元数据是否一致。`gap_filter` 与 `gap_bucket` 是独立开关；`set_use_gap` 只是同时设置两者的简写。

## SubseqHash2 适配边界

本包依据 [Shao-Group 官方 SubseqHash2 实现](https://github.com/Shao-Group/subseqhash2)的 ABCk orders、前向/反向动态规划和 pivot 合并思路，提供 **20 种标准氨基酸 + X 的标量 Rust 适配**。默认参数为 `window=25, subseq=11, d=5, repeats=11, table_seed=42, chunk_size=500`，未在 Chainmer 测试标签上调优。

它不是官方四字符核酸 SIMD 二进制，不能与官方实现的运行时间直接等同。包含 X 的窗口不产种子。前向和反向 DP 均有朴素子序列枚举对照测试，分块与整段处理有一致性测试。

## 哈希与兼容性

哈希以 `u128` 存储，但两个 64 位分量来自同一个 64 位滚动状态，因此**有效状态最多 64 bit**，不得宣称具有 128-bit 熵或碰撞保证；锚点组合只使用子哈希低 64 位，不会进一步减少原状态的信息量。当前索引格式为 `chainmer-index/4`，旧版索引必须重建。`chainmer_keys_legacy` 仅保留 v0.4 单通道、未屏蔽 X 的键输出用于历史核对，不参与新索引和默认检索。

FASTA 使用 `BufRead` 逐行解析，避免额外常驻完整文件字符串；当前 API 仍返回完整序列映射，并非全流式蛋白组索引。索引保存先写同目录临时文件并 `flush/sync_all`；Unix 直接原子替换，Windows 使用备份与失败回滚避免覆盖失败丢失旧索引。

177 条基准序列、258,610 个 v0.4 Chainmer 键已逐项核对，Rust 与 Python 原型不一致的序列数为 0；该兼容性不代表 v0.5.1 与旧版产生相同键集或排名。

## 修复前 v0.5.0 反馈子集 smoke test（历史结果）

使用反馈包四个 `tree_input.fa` 展示树子集实际运行新默认检索，排除 QUERY 自身后：古菌 7/7 参照全部进入候选且未报出其他候选；细菌 18/19；真菌 17/17 位于前 17；病毒 1/1。细菌弱参照 `WP_407854620.1` 仍未召回。**这些文件不是原始全库，结果只验证真实序列上的运行和排序行为，不是全库 AP、FP/FN 或性能结论。**

## 修复版基准摘要（177 序列复评）

同一 177 条 [NiFe] 氢化酶大亚基数据、留属标签和无向去重配对下，修复版 v0.5.1 的数据库 IDF 使用全库无标签统计，分类阈值仍仅由训练折确定：

| 方法 | Task3 AP（5 折） | Task4 AP（4 折） |
|---|---:|---:|
| Chainmer v0.5.1 双通道 | **0.825 ± 0.042** | **0.915 ± 0.056** |
| Chainmer v0.5.1 class 单通道 | 0.810 ± 0.034 | 0.903 ± 0.070 |
| Chainmer v0.5.1 raw 单通道 | 0.740 ± 0.041 | 0.863 ± 0.086 |

修复前双通道为 0.817 / 0.913，本轮复评为 0.825 / 0.915；该同数据对照不能拆解单项修复贡献，也不替代跨家族独立验证。双通道仍未超过 v0.4 SWPN 历史结果（0.916 / 0.949）或 Task4 的 BLASTp（0.990）。严格链化后的单次全部查询耗时由 1.923 s 增至 4.546 s。完整口径、修复前后表格、逐折结果与限制见 [BENCHMARK.md](BENCHMARK.md)。

## 许可

MIT，见 [LICENSE](LICENSE)。
