<p align="center">
  <img src="Chainmer_Logo.png" alt="Chainmer Logo" width="420">
</p>

# Chainmer v1.0.3

**语言 / Language：** 中文 | [English](README.md)

Chainmer 是一个用于**蛋白同源候选检索与预过滤**的 Rust 库和命令行工具。v1.0.0 将原来的"单通道类等价键计数"替换为 **raw/class 双通道位置种子、数据库 IDF 校准和共线链评分**；**v1.0.2 将默认间距桶宽 `delta` 由 15 调整为 5**（gap 分桶细化 3 倍，减少旁系偶然匹配），其余流水线与 v1.0.1 一致。旧计数与 SWPN 接口保留作消融对照。

> Seed Hit 只表示候选共享种子，不构成同源或直系同源判定。直系同源需要基因树—物种树协调等后续分析；HydDB 亚群标签也不等同于直系同源群。

## 核心组成

1. **双通道召回**：raw 通道使用原始氨基酸 `kc=5`，class 通道使用替换等价类 `class_kc=8`，两种键使用不同哈希盐；
2. **位置化锚点对**：保存查询与目标的双锚点位置，不再把种子仅视为无序键集合；
3. **数据库频率校准**：以当前数据库 document frequency 计算 IDF，高频键自动降权；
4. **重复贡献约束**：同一键在一个候选中的保留位置命中共同分摊一份 IDF 权重；每键以稳定哈希 bottom-k 最多保留 64 个代表命中；
5. **共线链评分**：在错开的对角线分桶内，以四个锚点端点在查询和目标中均严格递增为兼容条件求最优链；
6. **未知残基屏蔽**：X、B、Z、J、O、`*` 和其他非标准符号清洗为 X；覆盖 X 的 k-mer 不参与锚点选择；
7. **可复现索引**：记录实际参数、位置 postings、键频率、序列 IDF 权重、算法版本和数据库 BLAKE3 摘要。

## 安装与验证

需要 Rust 1.85 或更高版本。

```bash
cargo install --path . --locked
cargo test
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

## 快速开始

```bash
chainmer demo                                   # 合成序列演示，非生物学验证
chainmer seed -i db.faa -m chainmer             # 统计每条序列唯一种子键
chainmer seed -i db.faa -m chainmer --keys -o chainmer.tsv   # 导出 TSV（论文键口径）
chainmer index -d db.faa -o db.chainmer.idx --kc 5 --class-kc 8 --gmin 30 --gmax 300 --delta 5
chainmer search -d db.faa -q query.faa -x db.chainmer.idx --top 5   # 默认体制自适应判别
chainmer search -d db.faa -q query.faa --native --top 5             # 回退原生共线链打分
chainmer compare -i db.faa -o seed_compare.tsv  # 同二进制内种子数量/提取耗时比较
```

支持的 `seed --method`：`chainmer`、`kmer`、`spaced`、`minimizer`、`syncmer`、`strobemer`、`ckmer`、`subseqhash2`。syncmer 基线为 **Closed** 定义；spaced pattern 至少含一个 `1`；`simplified-strobemer` 为 fixed-offset anchor-pair（非标准 randstrobe/minstrobe）。

## 密度混合种子通道（v1.0.3 起为默认键集）

默认键集在 RAW5+CLASS8 锚点对键之上并入间距 k-mer 键：两掩码 `111101110111` / `111011010010111`，raw 与替换等价类双轨、步长 1、独立盐值命名空间。14 家族基因对真值（n=162）实测 Top-1 0.796→0.821、FN 0.199→0.024、MRR 0.861→0.890（七法第一），FP 0.41→0.46；对 kmer_diamond 呈「Top-1 差 1 查询、漏检 1/11」取舍。`ChainmerParams` 的 `dense_masks/dense_stride/dense_class` 可关闭该通道还原纯锚点对键集；索引 v7 内置 per-target 加密键表。

## 体制自适应判别（v1.0.3 起为 search 底层默认）

```bash
chainmer search -d homologs.faa -q queries.faa --top 5
# 可选旋钮：--lambda 1.0（背景扣除强度） --beta 1.0（子序列通道权重） --div-threshold 0.8
# --native 回退原生共线链打分；--weights 加权路径不受影响
chainmer adaptive -d homologs.faa -q queries.faa -o adaptive.tsv   # 全排序 TSV 导出（暴力参照路径）
```

链式种子召回层不变（接受门仍为共线键 IDF-Jaccard `Jc>0`），排序得分由体制门控：同源集间距直方图最小两两余弦 <0.8 判分歧体制（走原生 Jc）；保守体制叠加诊断子序列通道，高旁系风险者再叠加背景扣除（`score = Jc + β·Js − λ·Jbg`）。输入仅 FASTA（ID 约定 `species__accession`），不需要树/HMM/注释；判别层数据（链键表/子序列键表/架构直方图）内置于索引 v7，旧索引需重建。评分为 f64 精度；8 家族基因对真值实测与代价见 CHANGELOG v1.0.3（β-tubulin、ICD 存在真实回归，代价敏感场景用 `--native`）。`adaptive` 子命令保留为全排序导出（`query/regime/rank/id/score/gate`，FASTA 文件序平局裁决）；默认 search 的平局裁决按索引 catalog（ID 升序）。

## v1.0.2 默认参数

`kc=5`、`class_kc=8`、`w1=6`、间距窗 `[30,300]`、**`delta=5`**、每锚点最多 12 对。`delta` 由 v1.0.1 的 15 改为 5，是 30 配置三阶段扫描中唯一弱 Pareto 优于旧默认的配置（见 CHANGELOG v1.0.2）。修改 `kc`/`class_kc`/`delta` 会改变种子键与索引内容，**必须重建索引**。

## v1.0.2 基准（HydDB 2012 / 属水平留出 / KEGG 92，种子键检索口径）

口径：统一 IDF 加权 Jaccard 种子键集合评分（论文主口径），留一法／属水平留出；数据与全部方法输出见
[seed-baselines/data](https://github.com/rsj99-dev/Chainmer/tree/main/seed-baselines)（`data/hydrogenase_datasets/`、`data/results/`）。

| 任务 | Top-1（95% CI） | Macro-F1 | MRR | Recall@5 | AUC | 应答率 |
|---|---:|---:|---:|---:|---:|---:|
| HydDB 2012 全库 | **0.9945**（0.991–0.998） | 0.9896 | 0.9954 | 0.9955 | 0.9244 | 1.000 |
| 属水平留出 | **0.9275**（0.900–0.950） | 0.9013 | 0.9476 | 0.9675 | 0.8896 | 1.000 |
| KEGG 92 外部家族 | **1.0000** | 1.0000 | 1.0000 | 1.0000 | 0.8478 | 1.000 |

全量 2012 逐查询旁系分析：**假阳性（旁系多检）295（14.7%）、假阴性（直系漏检）0**；中位密度 **6.049 键/残基**、每序列中位键数 3086。

留出集 McNemar（Chainmer 对其余）：对 k-mer (BLAST+) 16:18（p=0.864，不显著）；对 Minimizer 43:6、Syncmer 41:17（p=0.0022）、SubseqHash2 53:11、Strobemer 174:2、k-mer (DIAMOND) 40:14（p=0.0005）均显著。即 Chainmer 显著超过除穷举 k-mer (BLAST+) 外的全部通用种子。

> 注意：上表为**种子键检索**口径（论文主对比）。本包 `search` 的共线链评分是另一口径；v1.0.0/v1.0.1 在 177 序列上的共线链 AP 历史结果见 [BENCHMARK.md](BENCHMARK.md) 历史附录与 CHANGELOG，不应与上表混读。

## Rust API

- `chainmer_seeds` / `chainmer_keys` / `ChainmerParams`：双通道位置种子与键生成；
- `SeedIndex`：位置索引、频率统计、保存、加载和数据库校验；
- `search` / `CandidateScore`：IDF 归一化共线链检索、分通道诊断分及位置命中截断统计；
- `search_count` / `search_with_keys`：旧式倒排与穷举计数对照；
- `search_weighted_index` / `search_weighted_with_keys`：加权检索；
- `subseqhash2_keys` / `SubseqHash2Params`：SubseqHash2 蛋白标量适配。

## SubseqHash2 适配边界

依据 [Shao-Group 官方实现](https://github.com/Shao-Group/subseqhash2)的 ABCk orders、前向/反向 DP 与 pivot 合并思路，提供 **20 种标准氨基酸 + X 的标量 Rust 适配**（默认 `window=25, subseq=11, d=5, repeats=11, table_seed=42, chunk_size=500`）。非官方四字符核酸 SIMD 二进制；含 X 窗口不产种子。

## 哈希与兼容性

哈希以 `u128` 存储但有效状态最多 64 bit，不得宣称 128-bit 熵。当前索引格式 `chainmer-index/7`（紧凑二进制 CSR＋判别层＋加密键表），旧索引必须重建。`chainmer_keys_legacy` 仅保留 v0.4 单通道键输出作历史核对。

## 许可

MIT，见 [LICENSE](LICENSE)。
