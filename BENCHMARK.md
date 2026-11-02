# Benchmark（可复现说明）

数据：177 条 [NiFe] 氢化酶大亚基 tips；留属五折；标签为初步树推断标签。
**数据不随代码包分发**；完整数据与脚本见研究仓库。

## 命令（研究仓库）

```bash
python tests/test_evoseed.py
python scripts/rev_evoseed_eval.py      # Task3/Task4 + 消融 + 索引成本
python scripts/rev_gap_ablation.py      # 间距语义消融
python scripts/rev_paper_figures.py
python scripts/make_docx_paper.py docs/EvoSeed_算法论文.md docs/EvoSeed_算法论文.docx
```

## Task 3（直系 vs 旁系；均值±SD）

| 方法 | AUPRC | AUROC | MCC | Recall@FPR10% |
|---|---|---|---|---|
| EvoSeed-SWPN | 0.911 ± 0.047 | 0.899 | 0.573 | 0.748 |
| simplified strobemer | 0.885 ± 0.077 | 0.910 | 0.219 | 0.609 |
| DIAMOND | 0.822 ± 0.045 | 0.828 | 0.405 | 0.308 |
| BLASTp | 0.812 ± 0.058 | 0.781 | 0.546 | 0.346 |

## Task 4（跨亚群硬负）

BLASTp/DIAMOND 0.995 > simplified strobemer 0.962 > EvoSeed-SWPN 0.945。

## 间距消融（Task 3，计数评分）

full 0.806 / filter_only 0.779 / no_gap 0.762 / bucket_only 0.786（AUPRC）。

## 说明

- `simplified strobemer` 是 fixed-offset anchor-pair 基线，不是标准 randstrobe/minstrobe。
- 负向键与 SWPN 权重仅由训练折拟合；MCC 阈值来自训练折。
- MMseqs2 因构建在结果导出阶段崩溃未纳入；未做内存剖析。
- `demo.py` 与 CLI demo 的序列为**合成示例**，不是生物学验证证据。
