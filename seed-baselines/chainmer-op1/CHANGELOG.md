# Changelog

## v1.0.0 — Chainmer-op1 首个独立优化发布

自 Chainmer v1.0.1（紧凑 CSR 索引 + 持久化缓存）分叉，class 通道五组重分组并完成 δ 曲线扫描优化。包名/二进制名 `chainmer-op1`，库名 `chainmer_op1`。

### 相对 Chainmer 的变更

1. **class 通道五组全覆盖**：`VIL→a`、`STQ→b`、`AGP→c`、`DEKRCH→d`、`FWYMN→e`；未知残基统一 f（6 类，字母表大小与原版相同、分组完全不同）。raw 通道、位置化锚点对、IDF 校准、共线链评分不变。
2. **默认间距桶 δ=3**：δ 曲线扫描（基线 δ=15，约束为 FN/FP 不劣化、holdout 不降）选定 δ=3。HydDB 全库 Top-1 0.9945→0.9955、FP 190→165、FN 1 不变；holdout 0.8775→0.8850；KEGG92 1.0。对原生 Chainmer（δ=5）的 holdout 0.9275 仍显著劣（McNemar 11:28，p=0.0095），原生综合最优的结论不变。
3. **缓存隔离**：默认缓存目录 `chainmer-op1`、环境变量 `CHAINMER_OP1_CACHE_DIR`、缓存身份标签 `Chainmer-op1-dual-chain-xmasked-VIL-STQ-AGP-DEKRCH-FWYMN`，与 Chainmer 及其他变体互不复用。
4. class 映射改变使种子键/DF/IDF 全部变化，**旧索引必须重建**；加载其他变体索引会因 class_map 校验失败而明确报错。回归测试按新映射重新锁定 golden 值。

### 继承自 Chainmer v1.0.1

- 索引格式 `chainmer-index/5`（紧凑二进制 CSR，含数据区 blake3 与尾部哨兵）；旧 v4/JSON 索引给出明确重建提示（`chainmer-op1 index`）。
- 持久化缓存：`search` 未传 `-x` 时按“格式+算法+版本+哈希方案+参数+数据库摘要”定位缓存，命中直接复用；`--rebuild-index`、`--cache-dir`、`--no-cache` 可用。
- Chainmer v1.0.0 的十二项修复（四端点严格递增链化、X 屏蔽、稳定哈希 bottom-k、索引完整性校验、原子写入等）全部继承。

### 验证

- `cargo build --release`、`cargo test --release` 全部测试通过（41 项）。
- 实测口径与数字见 [BENCHMARK.md](BENCHMARK.md)；基准数据不随包分发。