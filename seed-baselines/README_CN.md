# seed-baselines — 全部种子算法的 Rust 实现与对应基准数据

公开仓库：**https://github.com/rsj99-dev/Chainmer/tree/main/seed-baselines**

本目录包含论文《使用 Chainmer 种子算法帮助蛋白质找直系亲属》（Chainmer 种子基准）所用**全部八类种子算法**的
Rust 实现，以及**对应基准数据**：

| 算法 | crate / 模块 | 命令行 |
|---|---|---|
| k-mer (BLAST+) | 根 crate `src/seeds.rs::kmers` | `seed-baselines --method kmer --k 3 --stride 1` |
| k-mer (DIAMOND) | 根 crate `src/seeds.rs::spaced_kmers` | `seed-baselines --method spaced --masks 111101110111,111011010010111` |
| Minimizer | 根 crate `src/seeds.rs::minimizers` | `seed-baselines --method minimizer --k 7 --w 10` |
| Syncmer | 根 crate `src/seeds.rs::syncmers_with_mode`（闭合式） | `seed-baselines --method syncmer --k 7 --s 5` |
| Strobemer | 根 crate `src/seeds.rs::strobemers`（固定偏移） | `seed-baselines --method strobemer --k 7 --w1 9 --w2 32 --t 2` |
| SubseqHash2 | 根 crate `src/subseqhash2.rs`（蛋白标量） | `seed-baselines --method subseqhash2 --window 25 --subseq 11 --d 5 --repeats 11` |
| **Chainmer** | `chainmer/`（完整包 v1.0.2，δ=5） | `chainmer seed --method chainmer --keys …` |
| **Chainmer-op1** | `chainmer-op1/`（完整包 v1.0.0，五组全覆盖 class 通道，δ=3） | `chainmer-op1 seed --method chainmer --keys …` |

## 目录结构

```
seed-baselines/
├── Cargo.toml, src/, tests/      六基线 crate（lib seed_baselines + bin seed-baselines），
│                                 与论文导出逐字节兼容；6 项黄金标准回归测试
├── chainmer/                     Chainmer v1.0.2 完整包（双通道锚点对种子、索引/检索）
├── chainmer-op1/                 Chainmer-op1 v1.0.0 完整包（五组全覆盖 class 通道、δ=3）
├── data/
│   ├── hydrogenase_datasets/     基准输入数据：HydDB 2012、属水平留出、KEGG 92
│   │                             （136 文件，SHA-256 MANIFEST.tsv）
│   └── results/                  全部 8 方法（含 Chainmer 与 Chainmer-op1）的基准输出：
│                                 metrics_toolB/native、analysis、汇总 CSV、full2012 假阳性/假阴性、
│                                 McNemar、fig5 树案例、逐查询预测、gene_tree_rooted.nwk、reproduction_check.md
├── README.md / README_CN.md / CHANGELOG.md / LICENSE / .gitignore
```

## 为什么用 Rust、为什么关闭向量加速

- 跨平台兼容性：同一 Rust 代码库，跨操作系统/体系结构语义与编译行为一致，避免 C++ 未定义行为与编译器/模板分歧。
- 确定性计算：统一关闭 x86 CPU 向量加速指令集——纯标量可移植代码、不含 SIMD 内联指令（SubseqHash2 为蛋白标量适配），
  键集合与评分跨平台逐字节可复现。
- 工程可靠性：编译期借用与生命周期检查在编译阶段排除内存类错误。
- 实验一致性：Chainmer 本身即 Rust；全部基线共享同一哈希（splitmix64）、清洗（U→C、非标准→X）、去重与评分代码，
  把实现层误差降至最低。

## 数据说明

- `data/hydrogenase_datasets/` 为论文所用输入数据原件（SHA-256 校验）。
- `data/results/` 含**全部八方法（包括 Chainmer 与 Chainmer-op1）**的论文基准输出：准确率、宏平均 F1、平均倒数排名、
  前 5 位召回率、AUC、应答率、混淆矩阵、逐同源组召回、全库假阳性/假阴性、McNemar 检验、fig5 树案例、逐查询预测。
- 逐序列种子键 TSV 为 GB 级，**不**随仓库存储；可经上述命令行**逐字节**重导（已验证 57/57 与论文在用导出 SHA-256 全同，
  见 `data/results/reproduction_check.md`）。

## 构建与测试

```
cargo build --release            # 六基线 crate（根目录）
cargo test  --release            # 6 项黄金标准回归 + 单元测试
cd chainmer      && cargo build --release   # Chainmer v1.0.2
cd chainmer-op1  && cargo build --release   # Chainmer-op1 v1.0.0
```

中文路径下按各 crate 文档使用 rust-lld（link-self-contained） workaround。

## 许可

MIT（见 LICENSE）。
