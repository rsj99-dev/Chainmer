# EvoSeed（直系同源导向的蛋白种子）

EvoSeed 构建**面向直系同源候选检索**的蛋白种子，由三部分组成：

1. **六类氨基酸替换等价编码**：把 20 种氨基酸映射到 6 个物理化学类别，容忍保守替换；
2. **有序锚点对（anchor pair）**：在类串上取 minimizer 作为锚点，把两锚点的**相对间距分桶**写入哈希键；
3. **可选的正负对比特异性权重**：仅在训练折学习，`w(k) = ln((f_pos+1)/(f_neg+1))`。

> 重要：一次命中只是候选（Seed Hit）。直系同源**不能**由种子判定，必须由基因树—物种树协调确定。
> HydDB 亚群是分类参照，不等同于直系同源群。

## 目录结构

```
evoseed/
  seeds.py      替换类、哈希、k-mer/间隔/minimizer/syncmer/strobemer 工具
  core.py       evoseed_keys()：类锚点 -> 有序间距分桶锚点对键
  generator.py  generate()/generate_from_fasta()
  index.py      倒排索引构建/保存/读取
  search.py     候选检索（共享键计数）
  scoring.py    一致性得分、特异性权重、加权得分、负向键集
  extend.py     可选 BLOSUM62 复核打分
tests/test_evoseed.py
demo.py
```

## 快速开始

```bash
python demo.py
python tests/test_evoseed.py     # 18 项测试
```

核心仅依赖 Python 标准库；`biopython` 为可选依赖（仅 `score_pair_align` 使用）。

## 评分与负向排斥语义

- 计数评分：`search_with_keys()` / `candidate_rank()` 返回共享键数；
- SWPN 加权评分：`contrastive_weights()` 学习 `ln((f_pos+1)/(f_neg+1))`，
  `candidate_rank_weighted()`（或 `search_weighted()`）按权重和排序；
- 负向排斥：`max_reject` 为**允许命中的负向键数量上限**，命中数 > max_reject 才剔除；
  若传入 `negative_keys` 却未提供候选键集 `db_keys`，API 直接报错，不会静默失效。

## 已报告的关键结果

在由系统发育推导标签的 [NiFe] 氢化酶大亚基基准（177 条序列；留属五折交叉验证；权重仅在训练折学习）上：

- 直系 vs 旁系（Task 3）：EvoSeed-SWPN **AUPRC 0.911（95% CI 0.877–0.951）**，
  高于 strobemer（0.885）、DIAMOND（0.822）、k-mer（0.820）与 BLASTp（0.812）；
- **MCC 0.573**，高于间隔种子（0.556）与 BLASTp（0.546）；AUROC 0.899 略低于 strobemer（0.910）；
- 直系 vs 其他亚群硬负（Task 4）：BLASTp/DIAMOND（0.995）优于 EvoSeed-SWPN（0.945）；
- 索引成本较高：约 **1461 个键/序列**、171 个候选/查询（strobemer 为 103 与 13）；
- 同口径消融：无权 0.806 / 仅负向 0.793 / 泄漏 0.772 / **正负对比 0.911**。

完整表格、消融与局限见论文。

## 引用

若使用本代码，请引用随附论文（EvoSeed：面向直系同源候选检索的替换等价类有序锚点对种子）。

## 许可

MIT（见 LICENSE）。
