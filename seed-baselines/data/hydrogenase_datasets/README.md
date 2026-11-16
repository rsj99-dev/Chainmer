# [NiFe] 氢化酶大亚基测试数据集

整理日期：2026-11-09

来源：本工作区 `测试/` 与 `测试2/`（原目录保持原样，未移动、未删除任何文件）。

整理原则：

- 只收集**验证所需的输入**：序列、标签、划分、MSA、ref-only 参考模型；
- **不复制**索引（`.idx.json`）、训练模型（`eco_model*.json`）、运行结果（candidates/scores/eco/pair_accept 等）；
- 所有数据文件均为字节级复制；
- `02_family_msa/family_msas.tsv` 与 `03_genus_holdout/train_manifest.tsv` 因原始清单里的绝对路径存在编码乱码，已按新目录结构重新生成为相对路径；
- 顶层与子目录使用 ASCII 命名，规避本项目已知的非 ASCII 路径/编码问题；说明文字使用中文。

## 目录结构

```text
hydrogenase_datasets/
├── README.md
├── MANIFEST.tsv                     # 本目录全部文件的 SHA-256 与字节数
├── 01_small_validation/             # 小数据集验证（32 文件，约 1.7 MB）
│   ├── 01_three_family_representatives/
│   ├── 02_frhA_missing_rescue/
│   ├── 03_refonly_model/
│   ├── 04_kegg_external/
│   │   ├── id_lists/
│   │   └── kegg_variants/
│   ├── 05_negative_controls/
│   └── 06_search_reference_sets/
└── 02_large_validation/             # 大数据集验证（102 文件，约 12.7 MB）
    ├── 01_hyddb_2012/
    ├── 02_family_msa/
    ├── 03_genus_holdout/
    │   └── msa_train/
    └── 04_stress_and_subsets/
```

## 01_small_validation：小数据集验证

| 子目录 / 文件 | 序列数 | 内容与用途 |
|---|---:|---|
| `01_three_family_representatives/three_family_14.fasta` | 14 | 三家族核心测试集：6 条 MvhA、5 条 FrhA、3 条 VhtA；序列名带属名 |
| `01_three_family_representatives/MvhA.FASTA` | 1 | 单条 MvhA 代表序列（473 aa） |
| `01_three_family_representatives/FrhA.FASTA` | 1 | 单条 FrhA（437 aa） |
| `01_three_family_representatives/VhtA.FASTA` | 1 | 单条 VhtA（595 aa） |
| `01_three_family_representatives/FamillyIV.fasta` | 1 | 单条 Family IV（505 aa）；文件名沿用原始拼写 |
| `02_frhA_missing_rescue/missing_frh_real.faa` | 2 | 两条真实 FrhA，用于漏检补救（FN rescue）验证 |
| `03_refonly_model/` | 4 + 模型 | 四个家族 FASTA、`family_msas.tsv`、`refmodel.json`（约 1.4 MB，ref-only ECO 参考模型） |
| `04_kegg_external/kegg_MvhA.faa` | 27 | KEGG MvhA 家族外部验证集 |
| `04_kegg_external/kegg_FrhA.faa` | 60 | KEGG FrhA 家族外部验证集 |
| `04_kegg_external/kegg_VhtA.faa` | 5 | KEGG VhtA 家族外部验证集 |
| `04_kegg_external/kegg92.faa` + `kegg92_labels.tsv` + `kegg92_metadata.tsv` | 92 | KEGG 92 条外部验证集及家族标签/元数据 |
| `04_kegg_external/holdout_db.faa` / `holdout_query.faa` | 91 / 1 | KEGG 留出检索库与 query |
| `04_kegg_external/holdout2_db.faa` / `holdout2_query.faa` | 83 / 1 | 第二组 KEGG 留出检索库与 query |
| `04_kegg_external/id_lists/` | 4 个 TSV | `kegg_FrhA_ids`、`kegg_MvhA_ids`、`kegg_VhtA_ids`、`kegg_all_ids` |
| `04_kegg_external/kegg_variants/kegg_dedup.faa` | 84 | KEGG 去冗余集合 |
| `04_kegg_external/kegg_variants/kegg_dedup_query.faa` | 84 | 去冗余集合的 query 版本 |
| `04_kegg_external/kegg_variants/kegg_species_reduced.faa` | 60 | 按物种缩减后的集合 |
| `05_negative_controls/shuffled.faa` | 14 | 与三家族测试集一一对应的打乱序列，用作阴性对照 |
| `06_search_reference_sets/all18.faa` | 15 | `MvhA` 通用参考 + 14 条三家族序列 |
| `06_search_reference_sets/refs4.faa` | 1 | 历史 4 参考集合去重后的单条 MvhA |

## 02_large_validation：大数据集验证

| 子目录 / 文件 | 序列数 | 内容与用途 |
|---|---:|---|
| `01_hyddb_2012/hyddb.faa` | 2012 | HydDB [NiFe] 氢化酶大亚基全库 |
| `01_hyddb_2012/hyddb_labels.tsv` | 2012 | 每条序列的 group 标签，共 **29 个 group**（1a–1k、2a–2e、3a–3d、4a–4i） |
| `01_hyddb_2012/hyddb_metadata.tsv` | 2012 | group、物种、门、目、功能、氧敏感性、亚基数 |
| `01_hyddb_2012/hyddb_paralog_pairs.tsv` | — | 同物种内跨 group 的旁系对，用于旁系混淆分析 |
| `01_hyddb_2012/hyddb-results.csv` | 2012 | HydDB 原始导出表（分号分隔，含完整序列与注释），仅作溯源，不建议直接计算 |
| `02_family_msa/group_*.faa` | 29 组 | 每个 group 的原始序列 |
| `02_family_msa/group_*.msa.faa` | 29 组 | 每个 group 的 MSA |
| `02_family_msa/family_msas.tsv` | 29 行 | 相对路径清单：`1a<TAB>group_1a.msa.faa`（与 `hyddb_labels.tsv` 的 group 标签一致） |
| `03_genus_holdout/train.faa` | 1612 | 留属训练集 |
| `03_genus_holdout/test.faa` | 400 | 留属测试集 |
| `03_genus_holdout/test_labels.tsv` | 400 | 测试集 group 标签 |
| `03_genus_holdout/split.json` | — | 划分参数：`fraction=0.2`、`seed=42`、`n_train=1612`、`n_test=400` |
| `03_genus_holdout/msa_train/` | 29 组 | 训练集各 group 的 MSA |
| `03_genus_holdout/train_manifest.tsv` | 29 行 | 相对路径清单：`group<TAB>msa_train/x.msa.faa` |
| `04_stress_and_subsets/stress_db_3.faa` | 87 | 压力测试库（3 个物种规模） |
| `04_stress_and_subsets/stress_db_10.faa` | 279 | 压力测试库（10 个物种规模） |
| `04_stress_and_subsets/subset10.faa` | 10 | HydDB 10 条快速子集 |
| `04_stress_and_subsets/v15_query10.faa` | 10 | 外部 10 条 query 子集 |

## 未复制的内容（避免体积重复）

| 类型 | 原始位置示例 | 大小示例 | 说明 |
|---|---|---|---|
| 全库索引 | `测试2/data/hyddb.idx.json` | 约 830 MB | 派生文件，可用 `orthochainmer index --help` 重建 |
| 留属训练索引 | `测试2/holdout/train.idx.json` | 约 666 MB | 同上 |
| KEGG 索引 | `测试/kegg/*.idx.json` | 约 34 MB/个 | 同上 |
| ECO 模型 | `测试2/data/eco_model_v4_final.json` 等 | 约 9.6 MB/个 | 派生模型；本轮只保留了 ref-only 的 `refmodel.json` |
| 运行结果 | `测试2/data/v15_*.tsv`、`candidates/scores/pair_accept` 等 | — | 结果不是输入数据，保留在原目录 |
| 分析脚本 | `测试2/hyddb_*.py`、`测试/kegg_*.py` | — | 保留在原目录，便于复现 |

索引重建与参数请以 `orthochainmer index --help` 为准；本目录不包含任何索引缓存。

## 使用建议

- **小数据验证**：`three_family_14.fasta` 适合快速冒烟与家族分类回归；`shuffled.faa` 作阴性对照；`refonly_model/` 用于 ref-only 模型测试；KEGG 文件用于外部家族验证。
- **大数据验证**：`hyddb.faa` + `hyddb_labels.tsv` 做全库检索/评估；`03_genus_holdout/` 做留属泛化；`02_family_msa/` 用于 profile/MSA 相关流程；`04_stress_and_subsets/` 用于规模与耗时测试。
- 相关 CLI 入口：`orthochainmer search --help`、`orthochainmer profile build --help`、`orthochainmer orthology benchmark --help`。

## 数据语义与限制

1. HydDB 的 group 标签是数据库分组，用作 ortholog-group 代理，**不是经过独立实验验证的绝对真值**；不能据此宣称“零假阳性”。
2. `hyddb_paralog_pairs.tsv` 是同一物种内跨 group 的旁系候选对，不代表所有同物种拷贝都是旁系。
3. KEGG 标签是家族级（MvhA/FrhA/VhtA/Family IV），不是 pair 级 orthology 标签。
4. 所有序列按原始内容复制，未做去低复杂度、去 X、去冗余等处理；如需清洗请在副本上进行。
5. 原始目录 `测试/`、`测试2/` 仍是唯一权威来源；本目录是整理后的验证副本。

## 完整性

`MANIFEST.tsv` 列出了本目录下除自身以外所有文件的相对路径、SHA-256 与字节数，可用于校验复制完整性。
