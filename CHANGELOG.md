# Changelog

## v1.1.0 — 两段式检索、u64 键与持久化缓存（含密度混合键集）

- **新增 `stage2` 两段式检索**：第一段用 k=3 精确词倒排召回 top-N 候选；第二段仅对候选现场生成 Chainmer 链键（u64）与 Js 子序列键，复用 `adaptive_core` 判别（`Jc + β·Js − λ·Jbg`，外加 k-mer 兜底分）。Chainmer 键不建全库倒排，主索引仍为 k-mer 键表（约 3.7 MB），命中预计算缓存后无需重建键表。
- **14 家族质量（n=156，树真值）：** 召回 0.739、假阳性 16.0%、漏检 0.6%（对照：k-mer 0.740/21.2%/0%；全库 Chainmer 0.738/18.6%/1.9%）。
- **性能（710 条库 / 156 查询，单机实测）：** 命中缓存时批查询 **1.19 s（4 线程）/ 3.08 s（1 线程）**、峰值内存 **147 MB**；同机 DIAMOND 2.2.8 为 0.38/0.59 s、211 MB，blastp 2.17 为 5.65/6.36 s、24 MB（两者为官方 SIMD 构建，Chainmer 为标量实现）。
- **u64 键与缓存 `S2C2`：** stage2 内部键截断为 64 位（碰撞期望 ~1e-7），内存 272→147 MB；缓存文件 `<db>.stage2cache`（magic `S2C2` + db 长度 + mtime±2s 校验），命中加载 0.08 s，输出与重建逐字节一致（MD5 相同）。
- **新增 `--threads`**（预计算与查询并行，默认 4）与实验方法 `chainmer-class`（仅 CLASS8 锚点对键，用于消融；质量 0.645/21.2%/18.7%，不建议单独使用）。
- **CLI 精简：** 移除 `index` / `search` / `adaptive` 三个子命令（旧的全库链键倒排路径），仅保留 `seed`、`stage2`、`compare`；库内 `index`/`search` 模块仅保留供内部测试使用，不再作为 CLI 暴露。
- **无向量加速：** 全仓库为标量实现；`-C target-cpu=native` 实测零收益（1.23 s vs 1.23 s），热点为二分/HashMap/posting 指针追逐。

### v1.1.0 基础：密度混合种子通道并入默认键集（索引 v7）

- **算法**：默认键集从「RAW5+CLASS8 锚点对键」扩为「锚点对键 ∪ 两掩码间距键」（掩码 `111101110111` / `111011010010111`，raw 与替换等价类双轨、步长 1、独立盐值命名空间）；`ChainmerParams` 新增 `dense_masks` / `dense_stride` / `dense_class`（置空掩码可还原 v1.0.x 键集做消融）。锚点对通道、共线链打分与 adaptive 召回门（Jc>0）定义不变。
- **14 家族基因对真值实测（n=162，单协议 IDF-Jaccard）：** Top-1 0.7963→**0.8210**、FN 0.1989→**0.0239**（kmer_diamond 0.275 的 1/11）、MRR 0.8607→**0.8902**（七法第一）、FP 0.41→0.46。对 kmer_diamond（0.827/0.275）呈「Top-1 相差 1 个查询、漏检 1/11」的取舍；McNemar 对 kmer_diamond/syncmer/subseqhash2 p=1，对 strobemer p<1e-10。单家族：GLK 1.000、MvhA 0.833、McrA 0.778（七法第一）、Mtd/Hmd/FwdA/MtrA FN 清零、ICD FN 0.439→0.152。
- **索引格式 `chainmer-index/6` → `chainmer-index/7`：** 判别层新增 per-target 间距加密键表（16 字节键，与 CSR 位置 posting 分离）；索引键量约 +66%（键总数），建索引耗时同步上升；v6 及更旧索引明确拒绝并提示重建，缓存身份含格式＋版本自动失效。
- **确定性修复：** adaptive 背景扣除「较不相似拷贝」选择在 jc 平局时改按 ID 裁决（原按遍历序），使索引路径与暴力 `adaptive` 路径逐 ID、逐分严格一致（MtaA 三查询 top-20 全同）。
- **回归：** 53 项测试全绿（v7 roundtrip、dense 表独立重算对齐、索引/暴力路径一致）；`chainmer_keys_legacy`（v0.4 golden 306 键）保持；旧 v1.0.x 键集可由 `dense_masks=[]` 复现。
- **口径说明：** 预分析阶段树为快速口径（LG+G4＋1000 UFBoot）；HdrA 使用 BLAST 强同源子集（441 全集存档）；MtaA/FwdA/AcsA/MtrA 含注释审编。论文定稿时如需 ModelFinder/全量同源集复跑，指数会小幅变动。

## v1.0.3 — 体制自适应判别（adaptive）为 search 底层默认（索引 v6，--native 回退原生）

- **新增 `adaptive` 模块与子命令，并下沉为 `search` 默认路径（`search_adaptive`）：** 链式种子召回层一字不改（接受门仍为共线键 IDF-Jaccard Jc>0，`min_shared` 可收紧），仅对已召回候选做判别式重排 `score = Jc + β·Js − λ·Jbg`。Js＝诊断子序列通道（span20/weight8 四掩码、步长 2、blake2b-64）；Jbg＝判别式背景扣除（背景＝每多拷贝物种中较不相似于查询的拷贝，无监督旁系代理）。不需要树/HMM/注释等额外输入，仅依赖库内序列与 species__accession 命名约定。新增 `--native` 回退原生共线链打分，`--lambda/--beta/--div-threshold` 旋钮；`--weights` 加权路径不变；`adaptive` 子命令保留为全排序导出（FASTA 文件序平局裁决）。
- **体制门控：** 同源集 RAW 锚点间距直方图最小两两余弦 <0.8 ⇒ 分歧体制（走原生 Jc，保护 MvhA 类架构分歧甜区）；保守体制叠加 Js，再按旁系代理/top1 相似度比 >0.7 决定是否叠加背景扣除（保护饱和家族）。召回门与排序解耦，背景扣除不牺牲远缘召回（FN 不变）。判别层与暴力路径共用 `adaptive_core`，数值完全一致；默认 search 平局裁决按索引 catalog（ID 升序）。
- **索引格式升级 `chainmer-index/5` → `chainmer-index/6`：** 在紧凑二进制 CSR 之上内置判别层数据——每 target 的升序链键表、诊断子序列键表（blake2b-64）、RAW 锚点间距直方图（55 桶 f64）；meta 新增 `div_threshold` 与 `fam_divergent`（建索引时判定）。全部随索引持久化，检索时无需额外输入文件。旧 v5/v4/JSON 索引明确拒绝并提示重建；缓存身份含格式＋版本，旧缓存自动失效。
- **实测（8 家族基因对树真值 n=77，纯 f64 IDF-Jaccard 口径）：** 默认 search 与 Python f64 参照 77/77 查询体制标签/完整排序/得分完全一致（容差 1e-6）；合并 Top-1 0.727→0.753、fp_event 0.247→0.208、FN 0.203 不变；McrA 0.500→0.750（FP 0.500→0.250）、RecA 0.700→0.750、PG 0.556→0.667；MvhA/MtrA/GLK 持平。索引路径与暴力 adaptive 路径逐 ID 得分相同；recall@5 0.402 vs 0.413 的微差仅来自平局裁决序（文件序 vs ID 序）在同分候选间的排列，非数值差异。**诚实代价：** β-tubulin 1.000→0.917、ICD 0.917→0.833（子序列通道在部分饱和/低信号家族引入噪声）；PG recall@5 0.275→0.195。
- **重要口径发现：** 旧基准管线 scipy float32 CSR 在小分值下有显著舍入伪影（相对误差可达数十％），早期原型报告的“全家族无回归、合并 0.779”部分来自该伪影；本版采用 f64 评分，小幅增益与“无回归”类宣称均经 f64 复核。
- **回归保障：** 原生 seed/链式键/原生 search 函数零改动（历史 golden 键数 306 钉住）；新增 8 项 adaptive 单元测试＋4 项集成回归（原生键不变、解耦召回门、分歧体制保原生、v6 判别层 roundtrip、索引/暴力逐 ID 对齐）；全部 53 项测试通过。
- **定位说明：** adaptive 成为产品默认是工程集成决策，定位为有代价的工程缓解而非普适修复；科学代价不变（β-tubulin、ICD 真实回归），对代价敏感的场景用 `--native` 回退。

## v1.0.2 — 默认参数调优（delta 15→5）

- **默认 `delta` 15→5（`ChainmerParams::default()` 与 CLI `--delta` 缺省值同步）：** gap 分桶 `bucket = d / delta` 细化 3 倍，键更精确地编码锚点间距，减少旁系偶然匹配。30 配置三阶段扫描（kc×class_kc 4×4、特异性旋钮 10、delta 曲线 4）中唯一弱 Pareto 优于旧默认的配置。
- **实测（HydDB 2012 全库 n=2012，种子键集口径）：** Top-1 0.9935→0.9945（McNemar 0:2 零回归）；FP 313→295；FN 0→0；属留出 Top-1 0.9275 持平（McNemar 4:4，p=1.0）；KEGG92 Top-1 1.0000 持平；中位密度 6.049 keys/残基不变。
- **历史对照测试钉住旧参数：** `python_prototype_hash_and_key_parity` 显式使用 `delta: 15`，golden 值不随新默认漂移。

## v1.0.1 — 内存整改 A+B（紧凑 CSR + 持久化缓存）

- **紧凑索引（Milestone B）：** `SeedIndex` 改为 CSR：`keys: Vec<u128>` / `offsets: Vec<u64>` / `document_freq: Vec<u32>` / `postings: Vec<PackedPosting>`（`target_id/p1/p2` 各 u32，共 12 字节）/ `targets: Vec<TargetRecord>`。删除每条 posting 的 `String` 目标 ID、`BTreeMap<u128, Vec<SeedPosting>>` 嵌套结构、`target_keys` 全量 BTreeSet、逐 key DF BTreeSet，以及 `validate()` 中的 `observed` 全量 BTreeSet。
- **索引格式升级：** `chainmer-index/4`（JSON）→ **`chainmer-index/5`**（紧凑二进制 CSR，含数据区 blake3 与尾部哨兵）。旧 v4/JSON 索引给出明确迁移提示（`chainmer index` 重建）。
- **持久化缓存（Milestone A）：** `search` 未传 `-x` 时默认按“格式+算法+版本+哈希方案+参数+数据库摘要”定位持久缓存，命中直接复用，不再每次重建全库种子；新增 `--rebuild-index`、`--cache-dir`、`--no-cache`（缺省缓存目录可用 `CHAINMER_CACHE_DIR` 覆盖）。
- **流式读写：** `save()` 用 `BufWriter` 顺序写出并边写边算 blake3；`load()` 用尾部 meta + 分配前 section 边界校验 + 缓冲顺序读取，避免裸 `File` 逐字段读取的系统调用风暴。
- **实测（同一 2,012 序列 hyddb，keys/posting 数不变）：** 索引 830 MB→171.5 MB，建索引峰值 RSS 3,442 MB→314 MB，搜索峰值 RSS 2,663 MB→184 MB；cache hit 0.58 s。cache miss/hit/`-x`/`--no-cache` 四路径结果一致。
- **未完成：** 外存分块 builder（C）、mmap 只读后端（D）、v4 只读导入（改为明确重建提示）。

## v1.0.0

本次发布集中修复 Chainmer 的种子定义、共线链几何、模糊残基处理、位置命中保留、索引完整性和文件写入安全性问题。

### Bug 修复

1. **修复 Syncmer 的 Open/Closed 定义与 CLI 标签不一致。** 新增显式 `SyncmerMode::{Closed, Open}`；Closed 现在正确定义为最小 s-mer 位于首端或末端，Open 使用显式固定位置。CLI 基线更正为 Closed，并增加独立参考定义测试。旧布尔接口保留用于历史兼容。

2. **修复共线链可能连接端点顺序不一致命中的问题。** 旧实现只约束双锚点中心递增，可能串联交叉、反向或共享同角色端点的命中。v1.0.0 要求 `q1/q2/t1/t2` 四个端点均严格递增，并使用桶内动态规划求最大权重链；新增交叉端点、桶边界和独立子集枚举对照测试。

3. **修复 X 和其他非标准残基可制造虚假共享种子的问题。** Raw/Class 主通道现在在类别转换前共用有效性掩码；任何覆盖 X、B、Z、J、O、`*` 或其他非标准残基的 k-mer 均不参与锚点选择。`SecPolicy::C` 仍将 U 归一为 C，`SecPolicy::X` 将 U 按未知残基处理。legacy 键接口保留历史语义。

4. **修复每键只保留最先 64 个位置命中的遍历顺序偏倚。** 改为基于完整位置坐标的确定性稳定哈希 bottom-k，每键仍最多保留 64 个代表命中，但不再系统性偏向 N 端或 posting 前部。新增截断前命中数、保留后命中数和高频键统计，并验证 posting 反序后结果一致。

5. **修复损坏索引的冗余字段可能未经完整校验即参与检索的问题。** 现在从 postings 重新汇总每条序列的 unique key，校验 posting ID、`SequenceStats.keys`、`idf_weight` 和 document frequency 的一致性。零种子序列仍被允许，不再错误要求每条统计序列都必须出现在 postings 中。

6. **修复索引坐标缺少序列长度检查的问题。** 数据库校验现在按实际清洗后的序列长度检查 `p1/p2`，拒绝越界或极端位置，避免损坏索引进入链化。

7. **修复索引构建期间完整中间种子集合重复驻留的问题。** `SeedIndex::build` 改为逐序列生成种子并立即写入 postings，当前序列处理完成后释放中间集合；构建结果与物化参考实现逐项一致。

8. **修复极端参数和位置可能触发整数溢出或 panic 的风险。** 链化中心与对角线改用 `i128`，位置命中计数使用 checked arithmetic，IDF 计算避免整数加法溢出；新增 `usize::MAX` 参数和损坏位置测试。

9. **修复全 0 spaced pattern 生成无序列信息重复键的问题。** pattern 现在必须至少包含一个 `1`，`spaced_kmers("ACDE", "0000", 0)` 会明确返回参数错误。

10. **修复哈希有效状态说明可能引起误解的问题。** 保留原哈希和 golden vector，不破坏 Python 原型兼容性；文档明确 128-bit 只是存储宽度，两个分量来自同一 64-bit 状态，有效状态最多 64 bit。

11. **降低 FASTA 读取时的额外内存副本。** `read_fasta` 改用 `BufRead` 逐行解析，不再同时常驻完整文件字符串和解析后的序列映射。当前 API 仍返回完整序列映射，尚不是全流式蛋白组索引。

12. **修复索引写入中断可能留下部分 JSON 或破坏旧索引的问题。** 保存时先写同目录临时文件并执行 `flush/sync_all`；Unix 使用原子替换，Windows 使用备份与失败回滚，并验证已有索引可以安全覆盖。

### 兼容性变化

- 当前索引格式为 **`chainmer-index/4`**，索引算法标识为 `Chainmer-dual-chain-xmasked`。
- **v1.0.0 之前的索引必须重新构建**，不兼容版本会被明确拒绝，不会静默误读。
- 哈希方案保持不变；X 屏蔽和链化修复会改变种子或分数，因此不能把 v0.5.0 结果直接当作 v1.0.0 结果。

### 验证与同数据复评

- `cargo fmt --check`、`cargo clippy --locked --all-targets -- -D warnings` 通过。
- 40 项测试全部通过，release 构建、CLI 重复运行确定性、旧索引拒绝和 20 文件发布包验证通过。
- 同一 177 条 [NiFe] 氢化酶序列复评中，双通道 Task3 AP 为 0.825（v0.5.0：0.817），Task4 AP 为 0.915（v0.5.0：0.913）。
- 严格链化增加计算成本：全部 177 查询的单次耗时为 4.546 s（v0.5.0：1.923 s），directed hits 均为 5,033。

### 已知限制

同数据复评不能替代跨家族独立标签验证，也不能拆解单项修复的独立贡献。完整逐项消融、蛋白组规模峰值 RSS、全流式索引和 Linux 原生验证仍未完成。Chainmer 输出同源候选及排序，不构成直系同源判定。
