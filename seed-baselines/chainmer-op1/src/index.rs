use crate::core::{chainmer_seeds, ChainmerParams, PositionedSeed};
use crate::seeds::{class_map, clean_protein, SecPolicy};
use crate::{HASH_SCHEME, INDEX_FORMAT, VERSION};
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

/// v5 二进制魔数：`CHMNIDX` + 格式字节 `0x05`。
const FORMAT_MAGIC: [u8; 8] = *b"CHMNIDX\x05";
/// v5 尾部哨兵：`CHMNEND` + 格式字节 `0x05`。
const TRAILER_MAGIC: [u8; 8] = *b"CHMNEND\x05";
/// 头部固定长度：magic(8) + n_postings(8) + n_keys(8) + n_targets(8)。
const HEADER_LEN: u64 = 8 + 8 + 8 + 8;
/// 单条 packing posting 字节宽度：target_id(4) + p1(4) + p2(4)，小端。
const POSTING_BYTES: u64 = 12;
/// 尾部固定长度：meta_offset(8) + meta_len(8) + trailer_magic(8)。
const TRAILER_LEN: u64 = 24;
const MAX_META_BYTES: u64 = 256 * 1024 * 1024;
/// 单条 target 记录最小字节：id_len(4) + id(>=1) + seq_len(4) + unique(4) + idf(8)。
const MIN_TARGET_RECORD_BYTES: u64 = 4 + 1 + 4 + 4 + 8;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct DatabaseMeta {
    pub n_sequences: usize,
    pub content_blake3: String,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct IndexMeta {
    pub format: String,
    pub algorithm: String,
    pub version: String,
    pub params: ChainmerParams,
    pub class_map: BTreeMap<String, String>,
    pub hash_scheme: String,
    pub db: DatabaseMeta,
}

/// 紧凑位置 posting：target 用整数 catalog 下标而非字符串，坐标为清洗后序列坐标。
///
/// 字段宽度经 `const` 断言为 12 字节，与文件二进制布局一致；这消除了旧实现里
/// 每条 posting 各保存一份 `String` 目标 ID 的开销。
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct PackedPosting {
    pub target_id: u32,
    pub p1: u32,
    pub p2: u32,
}

const _: () = assert!(std::mem::size_of::<PackedPosting>() == POSTING_BYTES as usize);

/// 每条目标序列在 catalog 中只保存一次。
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct TargetRecord {
    pub id: String,
    pub sequence_len: u32,
    pub unique_key_count: u32,
    pub idf_weight: f64,
}

/// 紧凑 CSR 索引。
///
/// 键 `keys` 严格升序；slot 的 posting 范围是
/// `postings[offsets[slot]..offsets[slot + 1]]`，且按 `(target_id, p1, p2)` 升序、
/// 无重复。`document_freq[slot]` 与 `keys[slot]` 一一对应。`targets` 为按整数
/// `target_id` 索引的 catalog，等价于旧实现的 `sequence_stats`（多了 `sequence_len`）。
#[derive(Clone, Debug, PartialEq)]
pub struct SeedIndex {
    pub meta: IndexMeta,
    pub keys: Vec<u128>,
    pub offsets: Vec<u64>,
    pub document_freq: Vec<u32>,
    pub postings: Vec<PackedPosting>,
    pub targets: Vec<TargetRecord>,
}

#[derive(Serialize, Deserialize)]
struct BinaryMeta {
    schema: String,
    data_len: u64,
    data_blake3: String,
    n_keys: u64,
    n_postings: u64,
    n_targets: u64,
    meta: IndexMeta,
}

static WRITE_COUNTER: AtomicU64 = AtomicU64::new(0);

fn write_sibling(path: &Path, kind: &str) -> Result<PathBuf> {
    let name = path.file_name().context("索引路径必须包含文件名")?;
    let counter = WRITE_COUNTER.fetch_add(1, Ordering::Relaxed);
    Ok(path.with_file_name(format!(
        ".{}.{}.{}.{}.tmp",
        name.to_string_lossy(),
        std::process::id(),
        counter,
        kind
    )))
}

pub fn key_hex(key: u128) -> String {
    format!("{key:032x}")
}

pub fn parse_key_hex(text: &str) -> Result<u128> {
    if text.len() != 32 || !text.bytes().all(|b| b.is_ascii_hexdigit()) {
        bail!("无效的 128-bit 十六进制键：{text}");
    }
    Ok(u128::from_str_radix(text, 16)?)
}

pub fn database_meta(sequences: &BTreeMap<String, String>) -> DatabaseMeta {
    let mut hasher = blake3::Hasher::new();
    for (id, sequence) in sequences {
        hasher.update(&(id.len() as u64).to_le_bytes());
        hasher.update(id.as_bytes());
        hasher.update(&(sequence.len() as u64).to_le_bytes());
        hasher.update(sequence.as_bytes());
    }
    DatabaseMeta {
        n_sequences: sequences.len(),
        content_blake3: hasher.finalize().to_hex().to_string(),
    }
}

pub fn database_seeds(
    sequences: &BTreeMap<String, String>,
    params: &ChainmerParams,
) -> Result<BTreeMap<String, BTreeSet<PositionedSeed>>> {
    sequences
        .iter()
        .map(|(id, seq)| Ok((id.clone(), chainmer_seeds(seq, params)?)))
        .collect()
}

pub fn database_keys(
    sequences: &BTreeMap<String, String>,
    params: &ChainmerParams,
) -> Result<BTreeMap<String, BTreeSet<u128>>> {
    Ok(database_seeds(sequences, params)?
        .into_iter()
        .map(|(id, seeds)| (id, seeds.into_iter().map(|seed| seed.key).collect()))
        .collect())
}

/// 在任何大型数组分配之前，用 checked arithmetic 校验各 section 计数与数据区自洽：
/// 固定 section 不越界、catalog 能容纳 `n_targets` 条最小记录、所有计数可安全转 `usize`。
fn validate_layout(meta_offset: u64, n_postings: u64, n_keys: u64, n_targets: u64) -> Result<()> {
    let posting_bytes = n_postings
        .checked_mul(POSTING_BYTES)
        .context("索引 posting 区字节数溢出")?;
    let key_bytes = n_keys.checked_mul(16).context("索引 key 区字节数溢出")?;
    let offset_bytes = n_keys
        .checked_add(1)
        .and_then(|count| count.checked_mul(8))
        .context("索引 offsets 区字节数溢出")?;
    let df_bytes = n_keys
        .checked_mul(4)
        .context("索引 document_freq 区字节数溢出")?;
    let fixed = HEADER_LEN
        .checked_add(posting_bytes)
        .and_then(|total| total.checked_add(key_bytes))
        .and_then(|total| total.checked_add(offset_bytes))
        .and_then(|total| total.checked_add(df_bytes))
        .context("索引固定 section 总长度溢出")?;
    if fixed > meta_offset {
        bail!(
            "索引固定 section 长度 {fixed} 超过数据区长度 {meta_offset}（section 重叠或计数损坏）"
        );
    }
    let catalog_bytes = meta_offset - fixed;
    let min_catalog = n_targets
        .checked_mul(MIN_TARGET_RECORD_BYTES)
        .context("索引 target catalog 最小字节数溢出")?;
    if min_catalog > catalog_bytes {
        bail!("索引 target 数量 {n_targets} 与 catalog 区长度 {catalog_bytes} 不一致（至少需 {min_catalog} 字节）");
    }
    let _ = usize::try_from(n_postings).context("posting 数量超出平台范围")?;
    let _ = usize::try_from(n_keys).context("key 数量超出平台范围")?;
    let _ = usize::try_from(n_targets).context("target 数量超出平台范围")?;
    let _ = usize::try_from(meta_offset).context("meta 偏移超出平台范围")?;
    Ok(())
}

/// 统计本轮写入的数据区字节并同步更新 blake3 的写入包装器。
struct HashingWriter<W: Write> {
    inner: W,
    hasher: blake3::Hasher,
    written: u64,
}

impl<W: Write> Write for HashingWriter<W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let written = self.inner.write(buf)?;
        self.hasher.update(&buf[..written]);
        self.written += written as u64;
        Ok(written)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

fn read_u32<R: Read>(reader: &mut R) -> Result<u32> {
    let mut bytes = [0u8; 4];
    reader.read_exact(&mut bytes).context("索引文件被截断")?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_u64<R: Read>(reader: &mut R) -> Result<u64> {
    let mut bytes = [0u8; 8];
    reader.read_exact(&mut bytes).context("索引文件被截断")?;
    Ok(u64::from_le_bytes(bytes))
}

fn read_f64<R: Read>(reader: &mut R) -> Result<f64> {
    let mut bytes = [0u8; 8];
    reader.read_exact(&mut bytes).context("索引文件被截断")?;
    Ok(f64::from_le_bytes(bytes))
}

impl SeedIndex {
    pub fn build(sequences: &BTreeMap<String, String>, params: ChainmerParams) -> Result<Self> {
        params.validate()?;
        if sequences.is_empty() {
            bail!("空数据库不能建立索引");
        }
        if sequences.len() > u32::MAX as usize {
            bail!(
                "序列数 {} 超过 u32::MAX，PackedPosting.target_id 无法表示",
                sequences.len()
            );
        }
        let policy = if params.sec_policy == "C" {
            SecPolicy::C
        } else {
            SecPolicy::X
        };
        let n = sequences.len();
        let mut records: Vec<(u128, u32, u32, u32)> = Vec::new();
        let mut targets: Vec<TargetRecord> = Vec::with_capacity(n);
        for (tid, (id, sequence)) in sequences.iter().enumerate() {
            // 清洗后长度用于 validate_db 的坐标上界；chainmer_seeds 内部同样按此策略清洗。
            let cleaned_len = clean_protein(sequence, policy).len();
            if cleaned_len > u32::MAX as usize {
                bail!("序列 {id} 清洗后长度 {cleaned_len} 超过 u32::MAX");
            }
            targets.push(TargetRecord {
                id: id.clone(),
                sequence_len: cleaned_len as u32,
                unique_key_count: 0,
                idf_weight: 0.0,
            });
            let tid = tid as u32;
            // 只保留当前序列的中间种子，避免完整 database_seeds 与索引同时驻留。
            let seeds = chainmer_seeds(sequence, &params)?;
            for seed in seeds {
                if seed.p1 > u32::MAX as usize || seed.p2 > u32::MAX as usize {
                    bail!("序列 {id} 的种子坐标超过 u32::MAX");
                }
                records.push((seed.key, tid, seed.p1 as u32, seed.p2 as u32));
            }
        }
        // 全局排序：key 升序；同 key 内 (target_id, p1, p2) 升序。
        records.sort_unstable();

        let mut keys: Vec<u128> = Vec::new();
        let mut offsets: Vec<u64> = Vec::with_capacity(n + 1);
        let mut document_freq: Vec<u32> = Vec::new();
        let mut postings: Vec<PackedPosting> = Vec::with_capacity(records.len());
        let mut idx = 0usize;
        while idx < records.len() {
            let key = records[idx].0;
            let start = idx;
            while idx < records.len() && records[idx].0 == key {
                idx += 1;
            }
            let slice = &records[start..idx];
            // DF = 该 key 覆盖的不同 target 数量（posting 已按 target_id 分组）。
            let mut df: u32 = 0;
            let mut last: Option<u32> = None;
            for (_, tid, _, _) in slice {
                if last != Some(*tid) {
                    df = df.checked_add(1).context("document frequency 溢出")?;
                    last = Some(*tid);
                }
            }
            let idf = ((n as f64 + 1.0) / (df as f64 + 1.0)).ln();
            // 每个 target 对该 key 只计一次 unique key 与一份 IDF 权重。
            let mut last: Option<u32> = None;
            for (_, tid, _, _) in slice {
                if last != Some(*tid) {
                    let target = &mut targets[*tid as usize];
                    target.unique_key_count = target
                        .unique_key_count
                        .checked_add(1)
                        .context("unique key 计数溢出")?;
                    target.idf_weight += idf;
                    last = Some(*tid);
                }
            }
            keys.push(key);
            offsets.push(postings.len() as u64);
            document_freq.push(df);
            for (_, tid, p1, p2) in slice {
                postings.push(PackedPosting {
                    target_id: *tid,
                    p1: *p1,
                    p2: *p2,
                });
            }
        }
        offsets.push(postings.len() as u64);

        Ok(Self {
            meta: IndexMeta {
                format: INDEX_FORMAT.into(),
                algorithm: "Chainmer-dual-chain-xmasked".into(),
                version: VERSION.into(),
                params,
                class_map: class_map(),
                hash_scheme: HASH_SCHEME.into(),
                db: database_meta(sequences),
            },
            keys,
            offsets,
            document_freq,
            postings,
            targets,
        })
    }

    /// 键在排序数组中的 slot；缺失返回 `None`（语义等同旧 `BTreeMap::get`）。
    pub fn key_slot(&self, key: &u128) -> Option<usize> {
        self.keys.binary_search(key).ok()
    }

    /// 某 key 的 posting 切片（按 `(target_id, p1, p2)` 升序）；缺失返回 `None`。
    pub fn postings_for(&self, key: &u128) -> Option<&[PackedPosting]> {
        let slot = self.key_slot(key)?;
        let start = self.offsets[slot] as usize;
        let end = self.offsets[slot + 1] as usize;
        Some(&self.postings[start..end])
    }

    /// 某 key 的 document frequency；缺失按 0 处理（与旧语义一致）。
    pub fn document_frequency(&self, key: &u128) -> u32 {
        self.key_slot(key)
            .map(|slot| self.document_freq[slot])
            .unwrap_or(0)
    }

    pub fn idf(&self, key: &u128) -> f64 {
        let df = self.document_frequency(key) as f64;
        ((self.meta.db.n_sequences as f64 + 1.0) / (df + 1.0)).ln()
    }

    /// 校验索引与数据库一致：内容摘要、catalog（ID 与清洗后长度）以及 posting 坐标上界。
    ///
    /// 不再构造按序列的临时长度 Map；catalog 已按 `target_id`（等于 BTreeMap ID 顺序）
    /// 保存 ID 与清洗后长度，逐条对比即可。
    pub fn validate_db(&self, sequences: &BTreeMap<String, String>) -> Result<()> {
        self.validate()?;
        let actual = database_meta(sequences);
        if self.meta.db != actual {
            bail!("数据库与索引不匹配：序列 ID、顺序无关内容摘要或序列数不同");
        }
        if sequences.len() != self.targets.len() {
            bail!("索引 target 数量与数据库序列数不一致");
        }
        let policy = if self.meta.params.sec_policy == "C" {
            SecPolicy::C
        } else {
            SecPolicy::X
        };
        for ((id, sequence), record) in sequences.iter().zip(self.targets.iter()) {
            if &record.id != id {
                bail!("索引 catalog ID 与数据库不一致：{} != {id}", record.id);
            }
            let cleaned = clean_protein(sequence, policy).len();
            if record.sequence_len as usize != cleaned {
                bail!("序列 {id} 的 catalog 长度与清洗后长度不一致");
            }
        }
        for posting in &self.postings {
            let record = &self.targets[posting.target_id as usize];
            if posting.p1 >= record.sequence_len || posting.p2 >= record.sequence_len {
                bail!(
                    "索引 posting 坐标超出清洗后序列长度：{} ({}, {}) >= {}",
                    record.id,
                    posting.p1,
                    posting.p2,
                    record.sequence_len
                );
            }
        }
        Ok(())
    }

    /// 校验 CSR 结构、posting 排序/坐标、DF 与按 key 相邻统计的 target 级聚合。
    ///
    /// 临时内存为 O(targets)，而不是旧实现的 O(targets × 每序列 key 数) 全量 BTreeSet。
    pub fn validate(&self) -> Result<()> {
        if self.meta.format != INDEX_FORMAT {
            bail!("不支持的索引格式：{}", self.meta.format);
        }
        if self.meta.algorithm != "Chainmer-dual-chain-xmasked" || self.meta.version != VERSION {
            bail!(
                "索引算法或版本不兼容：{}/{}",
                self.meta.algorithm,
                self.meta.version
            );
        }
        if self.meta.hash_scheme != HASH_SCHEME {
            bail!("哈希方案不兼容：{}", self.meta.hash_scheme);
        }
        self.meta.params.validate()?;
        if self.meta.class_map != class_map() {
            bail!("索引残基类别映射与当前版本不一致");
        }
        if self.meta.db.n_sequences == 0 || self.targets.len() != self.meta.db.n_sequences {
            bail!("索引序列统计数量与数据库元数据不一致");
        }
        if self.keys.len() != self.document_freq.len() || self.offsets.len() != self.keys.len() + 1
        {
            bail!("CSR keys/offsets/document_freq 长度不一致");
        }
        if self.offsets.first().copied() != Some(0)
            || self.offsets.last().copied() != Some(self.postings.len() as u64)
        {
            bail!("CSR offsets 端点与 posting 总数不一致");
        }
        if self.keys.windows(2).any(|pair| pair[0] >= pair[1]) {
            bail!("索引 keys 未严格升序");
        }
        if self.offsets.windows(2).any(|pair| pair[0] > pair[1]) {
            bail!("索引 offsets 非单调");
        }
        if self
            .targets
            .iter()
            .any(|record| !record.idf_weight.is_finite() || record.idf_weight < 0.0)
        {
            bail!("索引包含无效的序列 IDF 权重");
        }
        let mut seen_ids: BTreeSet<&str> = BTreeSet::new();
        for record in &self.targets {
            if !seen_ids.insert(record.id.as_str()) {
                bail!("索引 target catalog 包含重复 ID：{}", record.id);
            }
        }

        // 按 key 扫描累积，仅用 O(targets) 的临时数组；等价于旧 observed 全量集合的判定。
        let mut unique = vec![0u64; self.targets.len()];
        let mut weight = vec![0f64; self.targets.len()];
        for slot in 0..self.keys.len() {
            let start = self.offsets[slot] as usize;
            let end = self.offsets[slot + 1] as usize;
            if start >= end {
                bail!("索引包含空 postings");
            }
            let slice = &self.postings[start..end];
            let mut df: u32 = 0;
            let mut previous_tid: Option<u32> = None;
            let mut previous: Option<(u32, u32, u32)> = None;
            for posting in slice {
                if posting.target_id as usize >= self.targets.len() {
                    bail!("索引 posting target_id 越界");
                }
                if posting.p2 <= posting.p1 {
                    bail!("索引 postings 包含无效锚点位置");
                }
                let tuple = (posting.target_id, posting.p1, posting.p2);
                if let Some(prev) = previous {
                    if tuple <= prev {
                        bail!("索引 postings 未按 (target_id, p1, p2) 严格升序或含重复位置记录");
                    }
                }
                previous = Some(tuple);
                if previous_tid != Some(posting.target_id) {
                    df = df.checked_add(1).context("document frequency 溢出")?;
                    previous_tid = Some(posting.target_id);
                }
            }
            if df == 0 || df as usize > self.meta.db.n_sequences {
                bail!("索引键频率超出数据库序列数范围");
            }
            if self.document_freq[slot] != df {
                bail!("索引键频率与 postings 不一致");
            }
            let idf = ((self.meta.db.n_sequences as f64 + 1.0) / (df as f64 + 1.0)).ln();
            let mut last: Option<u32> = None;
            for posting in slice {
                if last != Some(posting.target_id) {
                    let target = posting.target_id as usize;
                    unique[target] = unique[target]
                        .checked_add(1)
                        .context("unique key 计数溢出")?;
                    weight[target] += idf;
                    last = Some(posting.target_id);
                }
            }
        }
        for (slot, record) in self.targets.iter().enumerate() {
            if record.unique_key_count as u64 != unique[slot] {
                bail!("序列 {} 的 unique key 统计与 postings 不一致", record.id);
            }
            let expected = weight[slot];
            let tolerance = 1e-12 * expected.abs().max(1.0);
            if (record.idf_weight - expected).abs() > tolerance {
                bail!(
                    "序列 {} 的 IDF 权重与 postings/document frequency 不一致",
                    record.id
                );
            }
        }
        Ok(())
    }

    /// 流式写出 v5 二进制索引：`BufWriter` 顺序写，边写边算数据区 blake3，
    /// 完成后原子替换正式文件。
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        self.validate()?;
        let path = path.as_ref();
        if path.exists() && !path.is_file() {
            bail!("索引输出路径不是普通文件：{}", path.display());
        }
        let temporary = write_sibling(path, "write")?;
        let write_result = (|| -> Result<()> {
            let file = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temporary)?;
            let mut writer = HashingWriter {
                inner: BufWriter::new(file),
                hasher: blake3::Hasher::new(),
                written: 0,
            };
            writer.write_all(&FORMAT_MAGIC)?;
            writer.write_all(&(self.postings.len() as u64).to_le_bytes())?;
            writer.write_all(&(self.keys.len() as u64).to_le_bytes())?;
            writer.write_all(&(self.targets.len() as u64).to_le_bytes())?;
            for posting in &self.postings {
                writer.write_all(&posting.target_id.to_le_bytes())?;
                writer.write_all(&posting.p1.to_le_bytes())?;
                writer.write_all(&posting.p2.to_le_bytes())?;
            }
            for key in &self.keys {
                writer.write_all(&key.to_le_bytes())?;
            }
            for offset in &self.offsets {
                writer.write_all(&offset.to_le_bytes())?;
            }
            for df in &self.document_freq {
                writer.write_all(&df.to_le_bytes())?;
            }
            for record in &self.targets {
                let id = record.id.as_bytes();
                let id_len = u32::try_from(id.len()).context("target id 过长")?;
                writer.write_all(&id_len.to_le_bytes())?;
                writer.write_all(id)?;
                writer.write_all(&record.sequence_len.to_le_bytes())?;
                writer.write_all(&record.unique_key_count.to_le_bytes())?;
                writer.write_all(&record.idf_weight.to_le_bytes())?;
            }
            let meta_offset = writer.written;
            let data_blake3 = writer.hasher.clone().finalize().to_hex().to_string();
            let binary_meta = BinaryMeta {
                schema: INDEX_FORMAT.into(),
                data_len: meta_offset,
                data_blake3,
                n_keys: self.keys.len() as u64,
                n_postings: self.postings.len() as u64,
                n_targets: self.targets.len() as u64,
                meta: self.meta.clone(),
            };
            let meta_bytes = serde_json::to_vec(&binary_meta)?;
            writer.write_all(&meta_bytes)?;
            writer.write_all(&meta_offset.to_le_bytes())?;
            writer.write_all(&(meta_bytes.len() as u64).to_le_bytes())?;
            writer.write_all(&TRAILER_MAGIC)?;
            writer.inner.flush()?;
            writer.inner.get_ref().sync_all()?;
            Ok(())
        })();
        if let Err(error) = write_result {
            let _ = fs::remove_file(&temporary);
            return Err(error).with_context(|| format!("无法写入索引：{}", path.display()));
        }

        #[cfg(not(windows))]
        {
            if let Err(error) = fs::rename(&temporary, path) {
                let _ = fs::remove_file(&temporary);
                return Err(error).with_context(|| format!("无法替换索引：{}", path.display()));
            }
        }
        #[cfg(windows)]
        {
            if path.exists() {
                // std::fs::rename 在 Windows 不覆盖目标；备份用于失败回滚，避免丢失旧索引。
                let backup = write_sibling(path, "backup")?;
                if let Err(error) = fs::rename(path, &backup) {
                    let _ = fs::remove_file(&temporary);
                    return Err(error)
                        .with_context(|| format!("无法备份旧索引：{}", path.display()));
                }
                if let Err(error) = fs::rename(&temporary, path) {
                    let rollback = fs::rename(&backup, path);
                    let _ = fs::remove_file(&temporary);
                    return match rollback {
                        Ok(()) => Err(error).with_context(|| {
                            format!("无法替换索引，已恢复旧文件：{}", path.display())
                        }),
                        Err(rollback_error) => Err(error).with_context(|| {
                            format!(
                                "无法替换索引；旧文件保留在 {}，回滚失败：{}",
                                backup.display(),
                                rollback_error
                            )
                        }),
                    };
                }
                fs::remove_file(&backup).with_context(|| {
                    format!("新索引已写入，但无法清理备份：{}", backup.display())
                })?;
            } else if let Err(error) = fs::rename(&temporary, path) {
                let _ = fs::remove_file(&temporary);
                return Err(error).with_context(|| format!("无法写入索引：{}", path.display()));
            }
        }
        Ok(())
    }

    /// 读取 v5 二进制索引：先读尾部 meta，再做分配前的 section 边界校验与全程 blake3 校验，
    /// 最后顺序读取各 section。旧 v4/JSON 索引给出明确迁移提示而非静默误读。
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let mut file =
            File::open(path).with_context(|| format!("无法读取索引文件：{}", path.display()))?;
        let len = file.metadata()?.len();
        if len < HEADER_LEN + TRAILER_LEN {
            bail!("索引文件过短，不可能是完整的 v5 索引");
        }
        let mut magic = [0u8; 8];
        file.read_exact(&mut magic)?;
        if magic != FORMAT_MAGIC {
            bail!(
                "索引不是 chainmer-index/5 二进制格式；旧 v4/JSON 索引需用 `chainmer-op1 index` 重建"
            );
        }
        let n_postings = read_u64(&mut file)?;
        let n_keys = read_u64(&mut file)?;
        let n_targets = read_u64(&mut file)?;
        file.seek(SeekFrom::End(-(TRAILER_LEN as i64)))?;
        let meta_offset = read_u64(&mut file)?;
        let meta_len = read_u64(&mut file)?;
        let mut trailer_magic = [0u8; 8];
        file.read_exact(&mut trailer_magic)?;
        if trailer_magic != TRAILER_MAGIC {
            bail!("索引尾部标记损坏：文件被截断或以错误格式写入");
        }
        if meta_offset < HEADER_LEN || meta_offset > len {
            bail!("索引 meta 偏移越界：{meta_offset} > {len}");
        }
        let expected = meta_offset
            .checked_add(meta_len)
            .and_then(|value| value.checked_add(TRAILER_LEN))
            .context("索引 meta 长度溢出")?;
        if expected != len {
            bail!("索引长度不一致：声明 {expected}，实际 {len}");
        }
        if meta_len == 0 || meta_len > MAX_META_BYTES {
            bail!("索引 meta 长度非法：{meta_len}");
        }
        // 分配前先校验 section 边界与计数可转换性（P1-5）。
        validate_layout(meta_offset, n_postings, n_keys, n_targets)?;

        file.seek(SeekFrom::Start(meta_offset))?;
        let mut meta_bytes =
            vec![0u8; usize::try_from(meta_len).context("meta 长度超出平台范围")?];
        file.read_exact(&mut meta_bytes)?;
        let binary: BinaryMeta =
            serde_json::from_slice(&meta_bytes).context("索引 meta JSON 解析失败")?;
        if binary.schema != INDEX_FORMAT {
            bail!(
                "不支持的二进制索引 schema：{}（期望 {INDEX_FORMAT}）",
                binary.schema
            );
        }
        if binary.data_len != meta_offset {
            bail!("索引数据区长度与 meta 记录不一致");
        }
        if binary.n_keys != n_keys
            || binary.n_postings != n_postings
            || binary.n_targets != n_targets
        {
            bail!("索引头部计数与 meta 记录不一致");
        }

        // 流式校验数据区 [0, meta_offset) 的 blake3，复用固定大小缓冲。
        file.seek(SeekFrom::Start(0))?;
        let mut hasher = blake3::Hasher::new();
        let mut buffer = vec![0u8; 1 << 20];
        let mut remaining = meta_offset;
        while remaining > 0 {
            let take = remaining.min(buffer.len() as u64) as usize;
            file.read_exact(&mut buffer[..take])?;
            hasher.update(&buffer[..take]);
            remaining -= take as u64;
        }
        let actual = hasher.finalize().to_hex().to_string();
        if actual != binary.data_blake3 {
            bail!("索引内容校验失败：blake3 摘要不匹配（文件损坏或被修改）");
        }

        file.seek(SeekFrom::Start(HEADER_LEN))?;
        // 分段顺序读取用 1 MiB 缓冲：避免逐字段 read_exact 对裸 File 触发海量系统调用。
        let mut reader = BufReader::with_capacity(1 << 20, file);
        let mut postings = Vec::with_capacity(n_postings as usize);
        for _ in 0..n_postings {
            let target_id = read_u32(&mut reader)?;
            let p1 = read_u32(&mut reader)?;
            let p2 = read_u32(&mut reader)?;
            postings.push(PackedPosting { target_id, p1, p2 });
        }
        let mut keys = Vec::with_capacity(n_keys as usize);
        for _ in 0..n_keys {
            let mut bytes = [0u8; 16];
            reader.read_exact(&mut bytes)?;
            keys.push(u128::from_le_bytes(bytes));
        }
        let mut offsets = Vec::with_capacity(n_keys as usize + 1);
        for _ in 0..=n_keys {
            offsets.push(read_u64(&mut reader)?);
        }
        let mut document_freq = Vec::with_capacity(n_keys as usize);
        for _ in 0..n_keys {
            document_freq.push(read_u32(&mut reader)?);
        }
        let mut targets = Vec::with_capacity(n_targets as usize);
        for _ in 0..n_targets {
            let id_len = read_u32(&mut reader)?;
            if id_len == 0 {
                bail!("索引 target id 长度为零");
            }
            let mut id_bytes = vec![0u8; id_len as usize];
            reader.read_exact(&mut id_bytes)?;
            let id = String::from_utf8(id_bytes).context("索引 target id 不是合法 UTF-8")?;
            let sequence_len = read_u32(&mut reader)?;
            let unique_key_count = read_u32(&mut reader)?;
            let idf_weight = read_f64(&mut reader)?;
            targets.push(TargetRecord {
                id,
                sequence_len,
                unique_key_count,
                idf_weight,
            });
        }
        let consumed = reader.stream_position()?;
        if consumed != meta_offset {
            bail!("索引 target catalog 长度与 data_len 不一致：{consumed} != {meta_offset}");
        }

        let index = Self {
            meta: binary.meta,
            keys,
            offsets,
            document_freq,
            postings,
            targets,
        };
        index.validate()?;
        Ok(index)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn db() -> BTreeMap<String, String> {
        BTreeMap::from([("a".into(), "ACDEFGHIKLMNPQRSTVWY".repeat(8))])
    }

    fn temp_path(tag: &str) -> PathBuf {
        let counter = WRITE_COUNTER.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir().join(format!(
            "chainmer-compact-{tag}-{}-{}.ochidx",
            std::process::id(),
            counter
        ))
    }

    fn target_of(index: &SeedIndex, id: &str) -> u32 {
        index
            .targets
            .iter()
            .position(|record| record.id == id)
            .expect("target id 存在") as u32
    }

    #[test]
    fn corrupted_index_fields_are_rejected() {
        // target_id 越界
        let mut bad = SeedIndex::build(&db(), ChainmerParams::default()).unwrap();
        bad.postings[0].target_id = bad.targets.len() as u32;
        assert!(bad.validate().is_err());

        // unique_key_count 不一致
        let mut bad = SeedIndex::build(&db(), ChainmerParams::default()).unwrap();
        bad.targets[0].unique_key_count += 1;
        assert!(bad.validate().is_err());

        // idf_weight 不一致
        let mut bad = SeedIndex::build(&db(), ChainmerParams::default()).unwrap();
        bad.targets[0].idf_weight += 1.0;
        assert!(bad.validate().is_err());

        // document_freq 与 postings 不一致
        let mut bad = SeedIndex::build(&db(), ChainmerParams::default()).unwrap();
        bad.document_freq[0] = u32::MAX;
        assert!(bad.validate().is_err());

        // keys/offsets/document_freq 长度不一致
        let mut bad = SeedIndex::build(&db(), ChainmerParams::default()).unwrap();
        bad.keys.pop();
        assert!(bad.validate().is_err());

        // catalog 含重复 ID
        let mut bad = SeedIndex::build(&db(), ChainmerParams::default()).unwrap();
        bad.targets.push(TargetRecord {
            id: "a".into(),
            sequence_len: 1,
            unique_key_count: 0,
            idf_weight: 0.0,
        });
        assert!(bad.validate().is_err());
    }

    #[test]
    fn database_validation_rejects_out_of_bounds_positions_and_allows_zero_seed_sequences() {
        let database = db();
        let mut index = SeedIndex::build(&database, ChainmerParams::default()).unwrap();
        // 把最后一条 posting 的 p2 改为 u32::MAX：仍保持排序与 p2>p1，但超出清洗后长度。
        let last = index.postings.len() - 1;
        index.postings[last].p2 = u32::MAX;
        assert!(index.validate_db(&database).is_err());

        let database = BTreeMap::from([
            ("empty-seeds".into(), "A".into()),
            ("normal".into(), "ACDEFGHIKLMNPQRSTVWY".repeat(8)),
        ]);
        let index = SeedIndex::build(&database, ChainmerParams::default()).unwrap();
        index.validate().unwrap();
        index.validate_db(&database).unwrap();
        let empty = target_of(&index, "empty-seeds");
        assert_eq!(index.targets[empty as usize].unique_key_count, 0);
        assert!(index
            .postings
            .iter()
            .all(|posting| posting.target_id != empty));
    }

    #[test]
    fn incremental_build_matches_materialized_seed_reference() {
        let database = BTreeMap::from([
            ("a".into(), "ACDEFGHIKLMNPQRSTVWY".repeat(8)),
            ("b".into(), "ACDEFGHIKLMNPQRSTVWY".repeat(6)),
        ]);
        let params = ChainmerParams::default();
        let reference = database_seeds(&database, &params).unwrap();
        let index = SeedIndex::build(&database, params).unwrap();
        for (id, seeds) in reference {
            let target_id = target_of(&index, &id);
            for seed in seeds {
                let postings = index.postings_for(&seed.key).unwrap();
                assert!(postings.contains(&PackedPosting {
                    target_id,
                    p1: seed.p1 as u32,
                    p2: seed.p2 as u32,
                }));
            }
        }
        index.validate().unwrap();
        index.validate_db(&database).unwrap();
    }

    #[test]
    fn save_overwrites_an_existing_index_and_cleans_temporary_files() {
        let path = temp_path("roundtrip");
        let first = SeedIndex::build(&db(), ChainmerParams::default()).unwrap();
        first.save(&path).unwrap();
        let second_db = BTreeMap::from([("a".into(), "ACDEFGHIKLMNPQRSTVWY".repeat(9))]);
        let second = SeedIndex::build(&second_db, ChainmerParams::default()).unwrap();
        second.save(&path).unwrap();
        let loaded = SeedIndex::load(&path).unwrap();
        assert_eq!(loaded, second);
        assert_ne!(loaded.meta.db.content_blake3, first.meta.db.content_blake3);
        fs::remove_file(&path).unwrap();
    }

    #[test]
    fn load_rejects_legacy_json_and_non_v5_bytes() {
        let path = temp_path("legacy");
        fs::write(&path, b"{\"meta\":{}}").unwrap();
        assert!(SeedIndex::load(&path).is_err());
        fs::remove_file(&path).unwrap();
    }

    #[test]
    fn old_index_format_is_rejected() {
        let mut index = SeedIndex::build(&db(), ChainmerParams::default()).unwrap();
        index.meta.format = "chainmer-index/3".into();
        assert!(index.validate().is_err());
    }

    #[test]
    fn metadata_database_digest_and_positional_stats_are_real() {
        let mut params = ChainmerParams::default();
        params.set_use_gap(false);
        let index = SeedIndex::build(&db(), params).unwrap();
        assert!(!index.meta.params.gap_filter && !index.meta.params.gap_bucket);
        assert!(!index.meta.class_map.is_empty());
        assert!(index.postings.iter().all(|posting| posting.p2 > posting.p1));
        let a = target_of(&index, "a") as usize;
        assert_eq!(index.targets[a].unique_key_count as usize, index.keys.len());
        assert!(index.targets[a].idf_weight.is_finite());
        assert_eq!(index.targets[a].idf_weight, 0.0);
        let mut changed = db();
        changed.insert("b".into(), "ACDE".into());
        assert!(index.validate_db(&changed).is_err());
        let mut changed = db();
        changed.get_mut("a").unwrap().push('A');
        assert!(index.validate_db(&changed).is_err());
    }
}
