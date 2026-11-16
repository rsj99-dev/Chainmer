<p align="center">
  <img src="Chainmer_Logo.png" alt="Chainmer Logo" width="420">
</p>

# Chainmer-op1 v1.0.0

**语言 / Language：** 中文 | [English](README.md)

Chainmer-op1 是一个用于**蛋白同源候选检索与预过滤**的 Rust 库和命令行工具，自 Chainmer 分叉，采用五组全覆盖 class 通道并经 δ 曲线扫描优化默认间距桶（δ=3）；继承 Chainmer 的 **raw/class 双通道位置种子、数据库 IDF 校准和共线链评分**，旧计数与 SWPN 接口保留作消融对照。

> Seed Hit 只表示候选共享种子，不构成同源或直系同源判定。直系同源需要基因树—物种树协调等后续分析；HydDB 亚群标签也不等同于直系同源群。

## chainmer-op1 发布说明

本包是自 Chainmer v1.0.1 分叉的独立优化发布，包名与二进制名为 `chainmer-op1`（库名 `chainmer_op1`）：

- **class 通道五组覆盖全部 20 种标准氨基酸**：`VIL→a`、`STQ→b`、`AGP→c`、`DEKRCH→d`、`FWYMN→e`；未知残基统一为 f，共 6 类（与原版字母表大小相同，但分组完全不同）。
- raw 通道、位置化锚点对、IDF 校准、共线链评分与 Chainmer 完全相同；库 API 名称（`chainmer_seeds` 等）保持不变。
- 索引缓存目录（`chainmer-op1`）、环境变量 `CHAINMER_OP1_CACHE_DIR` 与缓存身份均与 Chainmer 及其他变体隔离，互不复用缓存。
- class 映射改变后全部种子键、DF/IDF 随之改变，**旧索引必须重建**；加载 Chainmer 或其他变体索引会因 class_map 校验失败而明确报错。
- 默认 `kc=5`、`class_kc=8` 沿用原版操作点；间距桶宽度默认 **δ=3**，由 δ 曲线扫描相对 δ=15 基线在 FN/FP 不劣化、holdout 不降约束下选定（HydDB Top-1 0.9945→0.9955、FP 190→165、holdout 0.8775→0.8850）。新分组下召回—候选量权衡会变化，部署前应用独立标签重新评估。

下文示例中的 `chainmer` 命令在本包中均读作 `chainmer-op1`（`seed --method chainmer` 为方法名，保持不变）。

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
chainmer-op1 demo

# 统计每条序列的唯一种子键
chainmer-op1 seed -i db.faa -m chainmer

# 导出 TSV：序列 ID + 32 位十六进制键
chainmer-op1 seed -i db.faa -m subseqhash2 --keys -o subseqhash2.tsv

# 构建索引并使用索引中保存的参数检索
chainmer-op1 index -d db.faa -o db.chainmer.json --kc 5 --class-kc 8 --gmin 30 --gmax 300
chainmer-op1 search -d db.faa -q query.faa -x db.chainmer.json --top 5

# 同一 Rust 二进制内的种子数量与单次提取耗时比较
chainmer-op1 compare -i db.faa -o seed_compare.tsv
```

支持的 `seed --method`：`chainmer`、`kmer`、`spaced`、`minimizer`、`syncmer`、`strobemer`、`ckmer`、`subseqhash2`。传统基线的固定参数见 `src/cli.rs`；syncmer 基线使用最小 s-mer 位于首端或末端的 **Closed** 定义，spaced pattern 至少要包含一个 `1`；`simplified-strobemer` 是 fixed-offset anchor-pair，不是标准 randstrobe/minstrobe。

默认 `search` 输出新的共线链归一化分数。`--raw` 和 `--class-only` 分别执行单通道消融。加权检索的 `--weights` JSON 沿用旧式加权重叠评分，以 32 位十六进制键为字段名、有限浮点数为值；它不是新的链化模型。`min_shared` 必须至少为 1；零权重共享键仍形成中性候选。加载已有索引时会校验数据库完整内容、冗余统计及清洗后序列坐标，而不是只比较序列 ID。

## 默认 k-mer 参数与调参边界

默认 `kc=5`、`class_kc=8` 是当前实现的**折中操作点**，不是理论上或跨数据集已经证明的最优组合。研究仓库在同一 177 条 [NiFe] 氢化酶序列上完成了 raw `kc=4/5/6` 与 class `class_kc=6/7/8/9` 的 12 组全因子实验，每组独立构建索引并运行全部查询 3 次。结果显示，raw 5 在每个受测 class 长度下都取得三种 raw 长度中最高的两种候选精确率；默认 5/8 没有被其他组合在目标对保留率、两种精确率和候选量上同时支配。该结果支持保留当前默认值，但不构成跨蛋白家族的最优性证明。

固定 raw `kc=5` 时，class 长度呈现明确的召回—候选量权衡：

| class_kc | 目标对保留率（%） | P_all（%） | P_Task3（%） | 候选/查询 |
|---|---:|---:|---:|---:|
| 6 | 88.48 | 15.84 | 71.97 | 69.073 |
| 7 | 80.35 | 22.67 | 74.11 | 43.819 |
| **8（默认）** | **74.31** | **33.48** | **80.98** | **27.435** |
| 9 | 67.55 | 39.46 | 83.69 | 21.164 |

目标对保留率基于 1,094 个无向 ortholog 对。P_all 以全部唯一候选对为分母，未标注对仍在分母中，因此是操作性候选精确率；P_Task3 仅使用具有互斥 ortholog/paralog 标签的候选对。候选量排除查询自身。更重视召回时可评估 class 7，更重视精确率或候选压缩时可评估 class 9；这些只是当前数据上的调参方向，部署到新数据库时仍应使用独立标签重新评估，不能只根据单次运行时间选择参数。

修改 `kc` 或 `class_kc` 会改变种子键、DF/IDF、候选和索引内容，必须**重新构建索引**；`search` 使用索引中保存的参数，不能只修改检索命令：

```bash
chainmer-op1 index -d db.faa -o db.kc5.class7.chainmer.json --kc 5 --class-kc 7
chainmer-op1 search -d db.faa -q query.faa -x db.kc5.class7.chainmer.json --top 5
```

研究仓库中的复现入口为 `python OrthologySeed/scripts/rev_kc_sensitivity.py --repeats 3`，输出为 `kc_sensitivity.tsv`、`kc_sensitivity_runs.tsv` 和 `kc_sensitivity_metadata.json`。基准序列与标签不随独立 Rust 代码包分发；完整实验口径和 12 组结果见论文补充数据 5。

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

哈希以 `u128` 存储，但两个 64 位分量来自同一个 64 位滚动状态，因此**有效状态最多 64 bit**，不得宣称具有 128-bit 熵或碰撞保证；锚点组合只使用子哈希低 64 位，不会进一步减少原状态的信息量。当前索引格式为 `chainmer-index/5`（紧凑二进制 CSR，整数 target_id），旧版索引必须重建。`chainmer_keys_legacy` 仅保留 v0.4 单通道、未屏蔽 X 的键输出用于历史核对，不参与新索引和默认检索。

FASTA 使用 `BufRead` 逐行解析，避免额外常驻完整文件字符串；当前 API 仍返回完整序列映射，并非全流式蛋白组索引。索引保存先写同目录临时文件并 `flush/sync_all`；Unix 直接原子替换，Windows 使用备份与失败回滚避免覆盖失败丢失旧索引。

177 条基准序列、258,610 个 v0.4 Chainmer 键已逐项核对，Rust 与 Python 原型不一致的序列数为 0；该兼容性不代表 v1.0.0 与旧版产生相同键集或排名。

## 基准摘要（Chainmer-op1 v1.0.0）

HydDB 2012 全库（2,012 条 [NiFe] 氢化酶大亚基、29 个亚群标签）上，Chainmer-op1 v1.0.0（kc=5、class_kc=8、δ=3）Top-1 为 0.9955；FP（最佳同物种旁系候选得分不低于最佳同物种直系候选的查询）165 个（8.2%），FN（漏检）1 个。种子密度 6.051 keys/aa（每序列中位 3,073 键、总键数 6,136,359），单次提取 2,012 条耗时 2.07 s（≈1.03 ms/序列）。留属评估 Top-1 0.8850，仍显著低于原生 Chainmer 0.9275（McNemar 11:28，p=0.0095）；KEGG 外部 92 条 Top-1 1.0。协议 A/B 下 Chainmer-op1 均用默认参数，两协议数字完全相同。完整口径、δ 扫描与限制见 [BENCHMARK.md](BENCHMARK.md)。

## 许可

MIT，见 [LICENSE](LICENSE)。
