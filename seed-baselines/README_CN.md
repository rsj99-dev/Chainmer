# seed-baselines v1.0.0

论文《使用 Chainmer 种子算法帮助蛋白质找直系亲属》（Chainmer 种子基准）六类基线种子算法的
**标量、无 SIMD 可移植 Rust 实现**：k-mer (BLAST+)、k-mer (DIAMOND)、Minimizer、Syncmer、
Strobemer、SubseqHash2。与论文导出工具（`paper_seed_benchmark/tools/seed_dump.rs`）逐字节兼容。

## 为什么用 Rust、为什么关闭向量加速

- **跨平台兼容性**：同一 Rust 代码库在操作系统与体系结构间语义与编译行为一致，避免 C++
  未定义行为、编译器与模板差异引入的实现分歧。
- **确定性计算**：统一关闭 x86 CPU 向量加速指令集——纯标量可移植代码、不含 SIMD intrinsics
  （SubseqHash2 为官方 C++ ABCk orders/前向反向 DP/pivot 合并的蛋白标量适配），
  键集合与评分跨平台逐字节可复现。
- **工程可靠性**：编译期借用与生命周期检查在编译阶段排除数据竞争、悬垂指针等内存类错误，
  全库长时导出与评测不因内存类 BUG 中断。
- **实验一致性**：Chainmer 内核本身即 Rust；本包与 Chainmer/op1 共享同一哈希（splitmix64）、
  清洗（U→C、非标准→X）、去重与评分代码，把实现层误差降至最低。

## 算法与参数（论文 toolB 主口径）

| 方法 | CLI | 参数来源 |
|---|---|---|
| k-mer (BLAST+) | `--method kmer --k 3 --stride 1` | blastp word_size=3 精确词（不复现 T=11 邻域词/two-hit） |
| k-mer (DIAMOND) | `--method spaced --masks 111101110111,111011010010111` | DIAMOND v2 默认双 weight-10 间隔形状（iedera）；按形状序号混盐 i×0x9E3779B97F4A7C15 取并集；不复现 Murphy10 缩减/双索引/频率过滤 |
| Minimizer | `--method minimizer --k 7 --w 10` | minimap2 k=15,w=10 按 f=log2(4)/log2(20)≈0.463 折算 |
| Syncmer | `--method syncmer --k 7 --s 5` | Weighted Syncmer (15,10) closed 定义折算；未含 weighted 排序 |
| Strobemer | `--method strobemer --k 7 --w1 9 --w2 32 --t 2` | StrobeMap 折算；fixed-offset 简化（非 randstrobe/minstrobe） |
| SubseqHash2 | `--method subseqhash2 --window 25 --subseq 11 --d 5 --repeats 11` | 蛋白标量适配；未用官方 SIMD |

## 构建与使用

```
cargo build --release
target/release/seed-baselines --input in.faa --output keys.tsv --method minimizer --k 7 --w 10
```

输出：`sequence_id<TAB>key_hex32<TAB>pos`（pos 为首次出现位置；集合型种子记 -1）。
逐序列按键去重保留首次位置；统计写 stderr（`method= seqs= residues= keys= elapsed_ms=`）。

## 测试

```
cargo test --release
```

六项回归 golden（`tests/fixtures/mini_*.tsv`，由论文在用导出二进制生成）与本包 CLI 输出逐字节比对，
另含清洗与参数校验单元测试。HydDB 2012 全库导出与论文在用 `seed_dump` 二进制逐字节一致
（见论文包 `results/reproduction_check.md`）。

## 许可

MIT（见 LICENSE）。
