# 复现一致性校验记录（reproduction_check）

强制全量重跑（run_all.py --force：57 键重导 → native/toolB 评测 → McNemar → 汇总 →
full2012 ×2 → fig5 缓存重算 → 中英 PNG/SVG → md/docx → 双 check）后与重跑前快照比对。

## 门禁项

- PASS metrics_native 键集合清洁（18）（new=18）
- PASS metrics_native 全字段一致 tol1e-9（18 条目）
- PASS metrics_toolB 键集合清洁（24）（new=24）
- PASS metrics_toolB 全字段一致 tol1e-9（24 条目）
- PASS preds 确定性：4 代表文件从键重算逐字节一致
- PASS preds 42 文件齐备（Top-1 经 verify_expected 重算≡旧 metrics）
- PASS full2012_seed_summary.json FP/FN 计数与比率一致（7 方法）
- PASS full2012_seed_summary_native.json FP/FN 计数与比率一致（7 方法）
- PASS full2012_seed_summary_toolB.json FP/FN 计数与比率一致（8 方法）
- PASS mcnemar_chainmer.json 一致
- PASS analysis summary native/toolB 行一致（新字段集⊆旧）
- PASS analysis mcnemar.toolB 一致
- PASS summary_native.csv 核心列一致（18 行）
- PASS summary_toolB.csv 核心列一致（24 行）
- PASS fig5_family_cases_toolB.json 一致
- PASS fig5_tree_cache_nine_toolB.npz 数组逐元素一致（32 数组）
- PASS group_case_cache_nine_toolB.npz 数组逐元素一致（16 数组）
- PASS 57 键 SHA-256 新旧全同
- PASS zh fig1_principle.png 与快照逐字节一致
- PASS zh fig3_main.png 与快照逐字节一致
- PASS zh fig4_errors.png 与快照逐字节一致
- PASS zh fig5_paralog.png 与快照逐字节一致
- PASS zh supp1_ranking.png 与快照逐字节一致
- PASS zh supp2_pipeline.png 与快照逐字节一致

## op1 导出 parity（前置门禁）

github/chainmer-op1 发布包新导出 ＝ 旧 chainmer3 变体 crate 导出 ＝ 论文在用键：
SHA-256 三方全同 `191E748AD1562C15D21074D17A21EC17232D37D84ED34716F1843EFC0554A168`
（hyddb toolB，294,419,346 字节）。

## 提取耗时新旧对照（ms/序列，HydDB；计时自然漂移属预期）

| 方法 | 旧 | 新 | 新舍入(2位) |
|---|---|---|---|
| kmer | 0.100 | 0.119 | 0.12 |
| kmer_diamond | 0.298 | 0.228 | 0.23 |
| minimizer | 0.150 | 0.156 | 0.16 |
| syncmer | 0.230 | 0.248 | 0.25 |
| strobemer | 0.150 | 0.155 | 0.15 |
| subseqhash2 | 58.510 | 32.513 | 32.51 |
| chainmer | 1.070 | 1.186 | 1.19 |
| chainmer3 | 1.030 | 1.051 | 1.05 |
| chainmer3(native) | 0.000 | 0.986 | 0.99 |

## 中文 PNG 与快照比对

- fig1_principle.png: identical
- fig2_cost_tradeoff.png: changed
- fig3_main.png: identical
- fig4_errors.png: identical
- fig5_paralog.png: identical
- supp1_ranking.png: identical
- supp2_pipeline.png: identical

## 论文 md 差异行（应仅附录 A 与表2 ms/措辞同步）

差异行数：24
```
-默认配置中位密度约 6.05 keys/aa（约为 k-mer (BLAST+) 的 6 倍、Minimizer 的 34 倍），HydDB 全库 v5 索引约 172 MB、提取约 1.1 ms/序列，定位为高召回、高成本的候选检索层。与 Strobemer 的区别在于：Chainmer 枚举多个下游锚点并用间距桶而非精确距离，对 indel 更宽容但键空间更大；其键含两个位置与间距，天然在召回时即携带共线信息，而非仅集合相似度。
+默认配置中位密度约 6.05 keys/aa（约为 k-mer (BLAST+) 的 6 倍、Minimizer 的 34 倍），HydDB 全库 v5 索引约 172 MB、提取约 1.2 ms/序列，定位为高召回、高成本的候选检索层。与 Strobemer 的区别在于：Chainmer 枚举多个下游锚点并用间距桶而非精确距离，对 indel 更宽容但键空间更大；其键含两个位置与间距，天然在召回时即携带共线信息，而非仅集合相似度。
-| k-mer (BLAST+) | k=3, stride=1 | blastp 默认 word_size=3 精确词（不复现邻域词/two-hit） | 0.942 | 0.10 |
-| k-mer (DIAMOND) | 2×weight-10 间隔种子（mask 111101110111 / 111011010010111） | DIAMOND v2 默认间隔种子（iedera 设计；不复现 Murphy10 缩减/双索引） | 1.950 | 0.30 |
-| Minimizer | k=7, w=10 | minimap2 k=15, w=10 折算 | 0.178 | 0.15 |
-| Syncmer | closed, k=7, s=5 | Weighted Syncmer (15,10) closed 折算 | 0.657 | 0.23 |
+| k-mer (BLAST+) | k=3, stride=1 | blastp 默认 word_size=3 精确词（不复现邻域词/two-hit） | 0.942 | 0.12 |
+| k-mer (DIAMOND) | 2×weight-10 间隔种子（mask 111101110111 / 111011010010111） | DIAMOND v2 默认间隔种子（iedera 设计；不复现 Murphy10 缩减/双索引） | 1.950 | 0.23 |
+| Minimizer | k=7, w=10 | minimap2 k=15, w=10 折算 | 0.178 | 0.16 |
+| Syncmer | closed, k=7, s=5 | Weighted Syncmer (15,10) closed 折算 | 0.657 | 0.25 |
-| SubseqHash2 | n=25, k=11, d=5, r=11（默认） | 默认 | 5.661 | 58.51 |
-| Chainmer | kc=5, class_kc=8, w1=6, gap[30,300], δ=5, 12 对（默认） | 默认 | 6.049 | 1.07 |
-
-*蛋白等效参数口径：核酸工具默认按字母表信息量折算（k 与 strobemer 偏移×f，f=log2(4)/log2(20)≈0.463，窗口保留）；k-mer 两组直接取蛋白工具默认、不折算。SubseqHash2/Chainmer 用默认参数；Chainmer-op1 为补充参数模式（源自补充文件的多家族旁系噪音分析），不入本表，参数（δ=3、五组全覆盖）与权衡见 3.5。提取时间为 HydDB 全库（2012 条）单线程总耗时除以序列数。*
+| SubseqHash2 | n=25, k=11, d=5, r=11（默认） | 默认 | 5.661 | 32.51 |
+| Chainmer | kc=5, class_kc=8, w1=6, gap[30,300], δ=5, 12 对（默认） | 默认 | 6.049 | 1.19 |
+
+*蛋白等效参数口径：核酸工具默认按字母表信息量折算（k 与 strobemer 偏移×f，f=log2(4)/log2(20)≈0.463，窗口保留）；k-mer 两组直接取蛋白工具默认、不折算。SubseqHash2/Chainmer 用默认参数；Chainmer-op1 为补充参数模式（源自补充文件的多家族旁系噪音分析），不入本表，参数（δ=3、五组全覆盖）与权衡见 3.5。提取时间为 HydDB 全库（2012 条）单线程总耗时除以序列数；Chainmer/op1 经发布包 CLI 导出，计时含进程启动开销（manifest timer=process-wall）。*
-成本—精度权衡首先体现在键密度与工程成本。密度差异极大：Strobemer 0.134、Minimizer 0.178、Syncmer 0.657、k-mer (BLAST+) 0.942、k-mer (DIAMOND) 1.950、SubseqHash2 5.66、Chainmer 6.05、Chainmer-op1 6.05 keys/aa，每序列中位键数 68–3086（图2A，八方法）。以 500 aa 计 Chainmer 约 3000 键，为 k-mer (BLAST+)（k=3）的 6 倍、Minimizer 的 34 倍——这是其以键空间换召回上限的设计选择。提取成本上 k-mer (BLAST+)/Minimizer/Syncmer/Strobemer 均 ≤0.25 ms/序列，k-mer (DIAMOND) 约 0.30、Chainmer 约 1.1、Chainmer-op1 约 1.0，SubseqHash2 因前向/反向 DP 约 58.5，高两个数量级。Chainmer 键集的 v5 二进制索引 137/172 MB、构建 1.9/2.4 s；相比 JSON 索引（666/830 MB、19.8/24.4 s）缩小约 4.8 倍、构建快约 10 倍，且逐行输出不变，使高密度选择在工程上可行。时间未做并行优化，仅作量级比较。
+成本—精度权衡首先体现在键密度与工程成本。密度差异极大：Strobemer 0.134、Minimizer 0.178、Syncmer 0.657、k-mer (BLAST+) 0.942、k-mer (DIAMOND) 1.950、SubseqHash2 5.66、Chainmer 6.05、Chainmer-op1 6.05 keys/aa，每序列中位键数 68–3086（图2A，八方法）。以 500 aa 计 Chainmer 约 3000 键，为 k-mer (BLAST+)（k=3）的 6 倍、Minimizer 的 34 倍——这是其以键空间换召回上限的设计选择。提取成本上 k-mer (BLAST+)/Minimizer/Syncmer/Strobemer 均 ≤0.25 ms/序列，k-mer (DIAMOND) 约 0.23、Chainmer 约 1.19、Chainmer-op1 约 1.05，SubseqHash2 因前向/反向 DP 约 32.5，高两个数量级。Chainmer 键集的 v5 二进制索引 137/172 MB、构建 1.9/2.4 s；相比 JSON 索引（666/830 MB、19.8/24.4 s）缩小约 4.8 倍、构建快约 10 倍，且逐行输出不变，使高密度选择在工程上可行。时间未做并行优化，仅作量级比较。
-chainmer-op1 与原生 Chainmer 的机理差异仅两处，其余流水线（RAW k=5 与 CLASS k=8 双通道、w1=6 Minimizer 锚点、[30,300] 间距窗、每锚点 12 对、同一哈希与清洗规则）逐行相同，故密度与提取成本持平（6.05 对 6.05 keys/aa、约 1.0 对 1.1 ms/序列）——差异不在索引规模，而在“哪些替换与距离漂移被视为同源信号”。其一，CLASS 通道等价类重分组：原生（AVLIMC、STNPG、DEQ、KRH、FYW）改为五组生化全覆盖（VIL、STQ、AGP、DEKRCH、FWYMN），替换容忍几何随之改变：原生容忍 A↔V/I/L 而 op1 不容忍，op1 容忍 D↔K 电荷反转与 F↔M、W↔N 而原生不容忍；近缘任务由 RAW 通道主导，重分组主要改变旁系间偶然共享 CLASS 键的模式，远缘任务则因容忍几何变化而减少共享键。其二，间距桶宽 δ 由 5 改 3：桶号=d/δ 整除，d∈[30,300] 的桶数由 55 增至 91，同桶距离漂移容忍由 4 aa 收紧到 2 aa。更细的桶即更强的共线约束：旁系锚点对即使相似、距离稍异即跨桶失效，偶然匹配键减少（FP 295→165）；但真实直系的 indel 漂移也更易跨桶丢键（FN 0→1、留出 0.9275→0.8850）。两版因此位于同一键设计“容忍—判别”权衡轴的两端：原生以宽桶与传统理化分组换远缘泛化与 indel 容忍，op1 以细桶与重分组换同物种旁系判别力（HydDB Top-1 0.9955 居首、FP 近乎减半），代价是远缘留出显著劣化；低 FP 与高全库排序不保证远缘泛化，选型须按任务距离谱取舍（4.3）。
+chainmer-op1 与原生 Chainmer 的机理差异仅两处，其余流水线（RAW k=5 与 CLASS k=8 双通道、w1=6 Minimizer 锚点、[30,300] 间距窗、每锚点 12 对、同一哈希与清洗规则）逐行相同，故密度与提取成本持平（6.05 对 6.05 keys/aa、约 1.05 对 1.19 ms/序列）——差异不在索引规模，而在“哪些替换与距离漂移被视为同源信号”。其一，CLASS 通道等价类重分组：原生（AVLIMC、STNPG、DEQ、KRH、FYW）改为五组生化全覆盖（VIL、STQ、AGP、DEKRCH、FWYMN），替换容忍几何随之改变：原生容忍 A↔V/I/L 而 op1 不容忍，op1 容忍 D↔K 电荷反转与 F↔M、W↔N 而原生不容忍；近缘任务由 RAW 通道主导，重分组主要改变旁系间偶然共享 CLASS 键的模式，远缘任务则因容忍几何变化而减少共享键。其二，间距桶宽 δ 由 5 改 3：桶号=d/δ 整除，d∈[30,300] 的桶数由 55 增至 91，同桶距离漂移容忍由 4 aa 收紧到 2 aa。更细的桶即更强的共线约束：旁系锚点对即使相似、距离稍异即跨桶失效，偶然匹配键减少（FP 295→165）；但真实直系的 indel 漂移也更易跨桶丢键（FN 0→1、留出 0.9275→0.8850）。两版因此位于同一键设计“容忍—判别”权衡轴的两端：原生以宽桶与传统理化分组换远缘泛化与 indel 容忍，op1 以细桶与重分组换同物种旁系判别力（HydDB Top-1 0.9955 居首、FP 近乎减半），代价是远缘留出显著劣化；低 FP 与高全库排序不保证远缘泛化，选型须按任务距离谱取舍（4.3）。
-脚本位于 paper_seed_benchmark/scripts/：bench.py（数据结构、IDF 加权 Jaccard 与指标）、dump_all.py（种子导出，含 TOOL_B 蛋白等效参数与 k-mer 双组）、run_bench.py（种子层）、full2012.py（全量 2012 逐 query 旁系多检/直系漏检，--seed-only --scope toolB）、mcnemar_chainmer.py（Chainmer 对通用种子 McNemar）、analyze.py（汇总）、fig5_trees.py（树案例，--family --scope toolB）、make_figs.py、build_paper.py。种子导出工具 tools/seed_dump.rs（--masks 支持多形状间隔种子并集、按形状序号混盐）需复制到随附 Chainmer 代码库 examples/（包内 vendor/）编译。中间结果在 results/，图在 figures/。复现：ASCII 路径下 cargo build --release --example seed_dump，依次运行 dump_all.py、run_bench.py --scopes toolB --methods kmer kmer_diamond minimizer syncmer strobemer subseqhash2 chainmer chainmer3、analyze.py、full2012.py --seed-only --scope toolB、fig5_trees.py --family --scope toolB、make_figs.py、build_paper.py、checks/verify_expected.py。随机性仅来自 bootstrap（seed=2026）与划分（seed=42），其余为确定性键集合运算。
+脚本位于 paper_seed_benchmark/scripts/：bench.py（数据结构、IDF 加权 Jaccard 与指标）、dump_all.py（种子导出：native 附录口径 + TOOL_B 蛋白等效主口径与 k-mer 双组 + Chainmer-op1 发布包键）、run_bench.py（种子层，--fresh 重建）、full2012.py（全量 2012 逐 query 旁系多检/直系漏检，--scope toolB/native）、mcnemar_chainmer.py（native Chainmer 对通用种子 McNemar）、analyze.py（汇总）、fig5_trees.py（树案例，--family --scope toolB）、make_figs.py（--lang both 中英双语 × PNG/SVG）、build_paper.py、meta.py（HydDB 元数据）；一键复现 run_all.py（--force 全量重导）。七类基线种子导出工具 tools/seed_dump.rs（--masks 支持多形状间隔种子并集、按形状序号混盐）由 tools/build_tools.ps1 复制到包内 vendor/（Chainmer v1.0.2 内核）examples/ 于 ASCII 路径编译；Chainmer 与 Chainmer-op1（内部键 chainmer3）键均由 github 发布包 seed --keys 导出（build_tools.ps1 编译 chainmer v1.0.2、build_op1.ps1 编译 op1 v1.0.0；op1 五组全覆盖 class 通道、δ=3），与历史键逐字节一致。各种子算法的 Rust 实现文件与函数映射见 种子算法源代码/README.md。中间结果在 results/，中文图在 figures/、英文图在 figures/en/（各含 PNG 与 SVG）。分步复现：依次运行 dump_all.py、run_bench.py --scopes native --methods kmer minimizer syncmer strobemer subseqhash2 chainmer --fresh、run_bench.py --scopes toolB --methods kmer kmer_diamond minimizer syncmer strobemer subseqhash2 chainmer chainmer3 --fresh、mcnemar_chainmer.py、analyze.py、full2012.py --scope native、full2012.py --scope toolB、fig5_trees.py --family --scope toolB、make_figs.py --lang both、build_paper.py、checks/verify_expected.py 与 checks/verify_inputs.py。随机性仅来自 bootstrap（seed=2026）与划分（seed=42），其余为确定性键集合运算。
```

结论：全部通过

## 记录后修订（2026-11-16，一致性门禁之后）

- 视觉修订：图2/3/5 方法色按用户要求重排——基线蓝绿系（BLAST+ 深青 #315C6B、DIAMOND 中青 #538A95、
  Minimizer 浅青 #8CB9BB、Syncmer 淡青 #C9DEDA、Strobemer 灰蓝 #607D8B、SubseqHash2 浅灰蓝 #A5BCC7），
  Chainmer/op1 黄棕系（#C6B992 / #7A6D45）。数据与版式不变；28 图文件与 docx 已重建。
- 新增 seed-baselines v1.0.0 包（github/seed-baselines）：六类基线标量 Rust 实现（无 x86 SIMD），
  六项回归 golden 全过；全库 HydDB 2012 六方法导出与论文在用 seed_dump 逐字节一致（6/6 SHA-256 全同）。
- 论文 §2.3 新增实现语言论证（Rust vs C++：跨平台兼容、关闭 x86 向量加速、编译期内存安全、同源一致性）；
  附录 A 增 seed-baselines 发布说明。以上均不影响本记录各门禁结论。
- 全文学术化润色（2026-11-16，用户要求）：pair 级→基因对水平（pair-wise）、group 级/同组→同源组水平、
  keys/aa→键/残基、ms/序列→毫秒/序列、FP/FN→假阳性/假阴性（首现定义）、query→查询序列、
  树证据→系统发育树证据、口径→参数方案、无命中率→无应答率、closed/fixed-offset/two-hit/indel 等改为中文术语＋英文括注；
  摘要/引言/讨论/结论/附录/图注全文统一。所有数字、表格数据、图与章节结构不变。
