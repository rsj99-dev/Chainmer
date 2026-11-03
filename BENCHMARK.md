# Benchmark 与复现说明

> 本文件数值均为 v0.4 历史结果。v0.5 已替换为双通道位置种子和共线链评分，尚未在本基准上重新评估；以下结果不得作为 v0.5 性能证据。

## 数据与评估口径

- 数据：177 条 [NiFe] 氢化酶大亚基序列；标签来自初步系统发育推断，不是直接的直系同源真值。
- 划分：按属哈希分配五折；任何触及测试属的配对均不进入该折训练集。
- 配对：无向去重，每对只评估一次。
- Task3：直系候选对 vs 旁系对，5 折；Task4：直系候选对 vs 跨亚群硬负，实际有 4 个有效折。
- SW/SWPN 权重和分类阈值仅由训练折拟合。下表 `Recall@FPR10%` 是测试 ROC 上 FPR 不超过 10% 的工作点，不是训练阈值迁移后的召回。
- AUPRC 为按同分组的非插值平均精度（AP），AUROC 也按同分组计算；全部同分时 AP=prevalence、AUROC=0.5。实现已与 scikit-learn 做随机并列分数对照。
- 数据不随代码包分发。以下数值来自实际运行，不包含官方 SubseqHash2 SIMD 结果。

## 判别性能

### Task3（均值 ± 折间 SD）

| 方法 | AUPRC | AUROC | MCC（训练折阈值） | Recall@FPR10% |
|---|---:|---:|---:|---:|
| EvoSeed-SWPN | 0.916 ± 0.045 | 0.897 ± 0.050 | 0.573 ± 0.167 | 0.748 ± 0.096 |
| EvoSeed 计数 | 0.802 ± 0.059 | 0.728 ± 0.065 | 0.353 ± 0.128 | 0.407 ± 0.149 |
| SubseqHash2 蛋白标量适配 | 0.719 ± 0.043 | 0.656 ± 0.027 | 0.300 ± 0.043 | 0.281 ± 0.033 |
| simplified strobemer | 0.684 ± 0.048 | 0.598 ± 0.058 | 0.219 | 0.244 |
| DIAMOND | 0.770 ± 0.042 | 0.727 | 0.405 | 0.308 |
| BLASTp | 0.818 ± 0.057 | 0.780 | 0.546 | 0.343 |

### Task4（均值 ± 折间 SD）

| 方法 | AUPRC | AUROC | MCC（训练折阈值） | Recall@FPR10% |
|---|---:|---:|---:|---:|
| EvoSeed-SWPN | 0.949 ± 0.045 | 0.900 ± 0.045 | 0.482 ± 0.094 | 0.740 ± 0.108 |
| EvoSeed 计数 | 0.919 ± 0.050 | 0.760 ± 0.051 | 0.332 ± 0.137 | 0.532 ± 0.079 |
| SubseqHash2 蛋白标量适配 | 0.867 ± 0.087 | 0.705 ± 0.029 | 0.316 ± 0.108 | 0.440 ± 0.070 |
| BLASTp | 0.990 ± 0.010 | 0.983 ± 0.012 | 0.531 | 0.966 |
| DIAMOND | 0.931 ± 0.051 | 0.843 ± 0.013 | 0.443 | 0.685 |
| simplified strobemer | 0.837 ± 0.112 | 0.637 ± 0.033 | 0.186 ± 0.064 | 0.318 ± 0.096 |

EvoSeed-SWPN 在 Task3 最高，BLASTp 在 Task4 最高。SubseqHash2 适配在两个任务上均未超过 EvoSeed 计数或 EvoSeed-SWPN；其 AUPRC/AUROC、训练阈值 MCC 和实际 FPR≤10% 召回均不领先 EvoSeed-SWPN。该单数据集结果不支持普遍性能排名。

## 种子数量、候选与提取成本

`rust_seed_compare.tsv` 来自同一 Windows 主机、同一 Rust release 二进制、177 条序列的一次顺序提取；耗时不是多次重复测量，不包含官方 C++ SIMD 对照，也未测量内存。

| 方法 | 键/序列 | 候选/查询（前 40 个查询） | 177 序列提取耗时 |
|---|---:|---:|---:|
| EvoSeed | 1461.1 | 171.4 | 29.260 ms |
| SubseqHash2 蛋白标量适配 | 2776.3 | 18.5 | 5387.955 ms |
| simplified strobemer | 102.9 | 13.2 | 9.139 ms |

SubseqHash2 的候选数明显少于 EvoSeed，但唯一键更多，且当前标量 DP 提取更慢。候选数来自统一研究评估脚本，提取耗时来自 Rust `compare`，两列不表示端到端检索时间。

## SubseqHash2 参数与实现边界

默认参数：`window=25, subseq=11, d=5, repeats=11, table_seed=42, chunk_size=500`。这些参数固定在代码中，未使用 EvoSeed 测试标签调优。

实现依据 Shao-Group 官方 C++ 的 ABCk 动态规划与 pivot 合并思路，字母表改为 20 种标准氨基酸 + X，使用标量 Rust。X 类似官方 N，包含 X 的窗口不产种子。它不是官方四字符核酸 SIMD 实现，不能把这里的耗时解释为官方 SubseqHash2 的性能。

## 复现命令

在发布目录构建：

```bash
cargo build --release
cargo test
```

在研究仓库中生成键集与种子层比较（路径按实际检出位置调整）：

```bash
evoseed seed -i analysis/alignment/nife_global_input.faa -m subseqhash2 --keys -o benchmark/orthology_tasks/subseqhash2_keys.tsv
evoseed compare -i analysis/alignment/nife_global_input.faa -o benchmark/orthology_tasks/rust_seed_compare.tsv
```

`rev_evoseed_eval.py` 支持通过 `EVOSEED_SUBSEQHASH2_KEYS` 环境变量读取上述 TSV，并沿用原有折、样本、权重和阈值流程。PowerShell 示例：

```powershell
$env:EVOSEED_SUBSEQHASH2_KEYS="benchmark/orthology_tasks/subseqhash2_keys.tsv"
python scripts/rev_evoseed_eval.py
```

主要输出：`evoseed_eval_folds.tsv`、`evoseed_eval_summary.tsv`、`seed_layer_comparison.tsv`、`rust_seed_compare.tsv`。

## 验证与限制

- Rust/Python v0.4 EvoSeed 全键集核验：177 条序列、258,610 个键、不一致序列数 0；v0.5 仅通过 `evoseed_keys_legacy` 保留该历史对照。
- v0.5 Rust 测试覆盖双通道位置种子、IDF 高频降权、目标额外键长度偏倚、共线链与散乱命中对照、位置索引存取、旧计数/穷举一致性、零权重候选、`min_shared`、负向排斥、重复 FASTA ID、参数错误，以及 SubseqHash2 DP 对照。
- 哈希只是 128-bit 存储，最多 64-bit 有效状态；不得宣称 128-bit 熵。
- 训练阈值迁移在其他分析中显示实际测试 FPR 可明显高于 10%，因此不能把训练 FPR 直接当作测试 FPR。
- 未做多次耗时统计、内存剖析或官方 SIMD 运行比较；MMseqs2 旧构建曾在结果导出阶段崩溃，未纳入本次比较。
