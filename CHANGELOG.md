# Changelog

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
