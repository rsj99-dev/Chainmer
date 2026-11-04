use crate::core::{chainmer_seeds, ChainmerParams, PositionedSeed};
use crate::seeds::{class_map, clean_protein, SecPolicy};
use crate::{HASH_SCHEME, INDEX_FORMAT, VERSION};
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

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

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
pub struct SeedPosting {
    pub id: String,
    pub p1: usize,
    pub p2: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SequenceStats {
    pub keys: usize,
    pub idf_weight: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SeedIndex {
    pub meta: IndexMeta,
    pub postings: BTreeMap<u128, Vec<SeedPosting>>,
    pub document_freq: BTreeMap<u128, usize>,
    pub sequence_stats: BTreeMap<String, SequenceStats>,
}

#[derive(Serialize, Deserialize)]
struct IndexFile {
    meta: IndexMeta,
    index: BTreeMap<String, Vec<SeedPosting>>,
    document_freq: BTreeMap<String, usize>,
    sequence_stats: BTreeMap<String, SequenceStats>,
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

impl SeedIndex {
    pub fn build(sequences: &BTreeMap<String, String>, params: ChainmerParams) -> Result<Self> {
        params.validate()?;
        if sequences.is_empty() {
            bail!("空数据库不能建立索引");
        }
        let mut postings: BTreeMap<u128, Vec<SeedPosting>> = BTreeMap::new();
        let mut target_keys: BTreeMap<String, BTreeSet<u128>> = BTreeMap::new();
        for (id, sequence) in sequences {
            // 只保留当前序列的中间种子，避免完整 database_seeds 与索引同时驻留。
            let seeds = chainmer_seeds(sequence, &params)?;
            let keys = target_keys.entry(id.clone()).or_default();
            for seed in seeds {
                keys.insert(seed.key);
                postings.entry(seed.key).or_default().push(SeedPosting {
                    id: id.clone(),
                    p1: seed.p1,
                    p2: seed.p2,
                });
            }
        }
        let n = sequences.len();
        let mut document_freq = BTreeMap::new();
        for (key, hits) in &postings {
            let df = hits
                .iter()
                .map(|hit| &hit.id)
                .collect::<BTreeSet<_>>()
                .len();
            document_freq.insert(*key, df);
        }
        let idf = |key: &u128| {
            let df = document_freq.get(key).copied().unwrap_or(0);
            ((n as f64 + 1.0) / (df as f64 + 1.0)).ln()
        };
        let sequence_stats = target_keys
            .into_iter()
            .map(|(id, keys)| {
                let stat = SequenceStats {
                    keys: keys.len(),
                    idf_weight: keys.iter().map(&idf).sum(),
                };
                (id, stat)
            })
            .collect();
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
            postings,
            document_freq,
            sequence_stats,
        })
    }

    pub fn idf(&self, key: &u128) -> f64 {
        let df = self.document_freq.get(key).copied().unwrap_or(0);
        ((self.meta.db.n_sequences as f64 + 1.0) / (df as f64 + 1.0)).ln()
    }

    pub fn validate_db(&self, sequences: &BTreeMap<String, String>) -> Result<()> {
        self.validate()?;
        let actual = database_meta(sequences);
        if self.meta.db != actual {
            bail!("数据库与索引不匹配：序列 ID、顺序无关内容摘要或序列数不同");
        }
        let policy = if self.meta.params.sec_policy == "C" {
            SecPolicy::C
        } else {
            SecPolicy::X
        };
        let lengths: BTreeMap<&str, usize> = sequences
            .iter()
            .map(|(id, sequence)| (id.as_str(), clean_protein(sequence, policy).len()))
            .collect();
        for hit in self.postings.values().flatten() {
            let length = lengths
                .get(hit.id.as_str())
                .copied()
                .context("posting ID 不在数据库中")?;
            if hit.p1 >= length || hit.p2 >= length {
                bail!(
                    "索引 posting 坐标超出清洗后序列长度：{} ({}, {}) >= {}",
                    hit.id,
                    hit.p1,
                    hit.p2,
                    length
                );
            }
        }
        Ok(())
    }

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
        if self.meta.db.n_sequences == 0 || self.sequence_stats.len() != self.meta.db.n_sequences {
            bail!("索引序列统计数量与数据库元数据不一致");
        }
        if self
            .sequence_stats
            .values()
            .any(|stat| !stat.idf_weight.is_finite() || stat.idf_weight < 0.0)
        {
            bail!("索引包含无效的序列 IDF 权重");
        }
        // 零种子序列合法，因此 posting ID 集合可以是统计 ID 集合的真子集。
        let mut observed: BTreeMap<String, BTreeSet<u128>> = self
            .sequence_stats
            .keys()
            .map(|id| (id.clone(), BTreeSet::new()))
            .collect();
        for (key, hits) in &self.postings {
            if hits.is_empty() {
                bail!("索引包含空 postings");
            }
            let unique: BTreeSet<_> = hits.iter().collect();
            if unique.len() != hits.len() {
                bail!("索引 postings 包含重复位置记录");
            }
            if hits.iter().any(|hit| hit.p2 <= hit.p1) {
                bail!("索引 postings 包含无效锚点位置");
            }
            let df = hits
                .iter()
                .map(|hit| &hit.id)
                .collect::<BTreeSet<_>>()
                .len();
            if df == 0 || df > self.meta.db.n_sequences {
                bail!("索引键频率超出数据库序列数范围");
            }
            if self.document_freq.get(key) != Some(&df) {
                bail!("索引键频率与 postings 不一致");
            }
            for hit in hits {
                let keys = observed
                    .get_mut(&hit.id)
                    .context("索引 posting 包含未知序列 ID")?;
                keys.insert(*key);
            }
        }
        if self.document_freq.len() != self.postings.len() {
            bail!("索引键频率表与 postings 大小不一致");
        }
        for (id, stat) in &self.sequence_stats {
            let keys = &observed[id];
            if stat.keys != keys.len() {
                bail!("序列 {id} 的 unique key 统计与 postings 不一致");
            }
            let expected: f64 = keys.iter().map(|key| self.idf(key)).sum();
            let tolerance = 1e-12 * expected.abs().max(1.0);
            if (stat.idf_weight - expected).abs() > tolerance {
                bail!("序列 {id} 的 IDF 权重与 postings/document frequency 不一致");
            }
        }
        Ok(())
    }

    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        self.validate()?;
        let path = path.as_ref();
        if path.exists() && !path.is_file() {
            bail!("索引输出路径不是普通文件：{}", path.display());
        }
        let file = IndexFile {
            meta: self.meta.clone(),
            index: self
                .postings
                .iter()
                .map(|(key, hits)| {
                    let mut hits = hits.clone();
                    hits.sort();
                    (key_hex(*key), hits)
                })
                .collect(),
            document_freq: self
                .document_freq
                .iter()
                .map(|(key, df)| (key_hex(*key), *df))
                .collect(),
            sequence_stats: self.sequence_stats.clone(),
        };
        let text = serde_json::to_string_pretty(&file)?;
        let temporary = write_sibling(path, "write")?;
        let write_result = (|| -> Result<()> {
            let mut output = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temporary)?;
            output.write_all(text.as_bytes())?;
            output.flush()?;
            output.sync_all()?;
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

    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let text = fs::read_to_string(path.as_ref())?;
        let file: IndexFile =
            serde_json::from_str(&text).context("索引必须符合 chainmer-index/4；旧版索引需重建")?;
        let mut postings = BTreeMap::new();
        for (key, hits) in file.index {
            let parsed = parse_key_hex(&key)?;
            if postings.insert(parsed, hits).is_some() {
                bail!("索引包含重复解析键：{key}");
            }
        }
        let mut document_freq = BTreeMap::new();
        for (key, df) in file.document_freq {
            let parsed = parse_key_hex(&key)?;
            if document_freq.insert(parsed, df).is_some() {
                bail!("索引频率表包含重复解析键：{key}");
            }
        }
        let index = Self {
            meta: file.meta,
            postings,
            document_freq,
            sequence_stats: file.sequence_stats,
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

    #[test]
    fn corrupted_redundant_index_fields_are_rejected() {
        let mut unknown_id = SeedIndex::build(&db(), ChainmerParams::default()).unwrap();
        unknown_id.postings.values_mut().next().unwrap()[0].id = "missing".into();
        assert!(unknown_id.validate().is_err());

        let mut wrong_keys = SeedIndex::build(&db(), ChainmerParams::default()).unwrap();
        wrong_keys.sequence_stats.get_mut("a").unwrap().keys += 1;
        assert!(wrong_keys.validate().is_err());

        let mut wrong_weight = SeedIndex::build(&db(), ChainmerParams::default()).unwrap();
        wrong_weight.sequence_stats.get_mut("a").unwrap().idf_weight += 1.0;
        assert!(wrong_weight.validate().is_err());

        let mut wrong_df = SeedIndex::build(&db(), ChainmerParams::default()).unwrap();
        let key = *wrong_df.document_freq.keys().next().unwrap();
        wrong_df.document_freq.insert(key, usize::MAX);
        assert!(wrong_df.validate().is_err());

        let mut missing_key = SeedIndex::build(&db(), ChainmerParams::default()).unwrap();
        let key = *missing_key.postings.keys().next().unwrap();
        missing_key.postings.remove(&key);
        missing_key.document_freq.remove(&key);
        assert!(missing_key.validate().is_err());

        let mut extra_stat = SeedIndex::build(&db(), ChainmerParams::default()).unwrap();
        extra_stat.sequence_stats.insert(
            "extra".into(),
            SequenceStats {
                keys: 0,
                idf_weight: 0.0,
            },
        );
        assert!(extra_stat.validate().is_err());
    }

    #[test]
    fn database_validation_rejects_extreme_positions_and_allows_zero_seed_sequences() {
        let database = db();
        let mut index = SeedIndex::build(&database, ChainmerParams::default()).unwrap();
        index.postings.values_mut().next().unwrap()[0].p2 = usize::MAX;
        assert!(index.validate_db(&database).is_err());

        let database = BTreeMap::from([
            ("empty-seeds".into(), "A".into()),
            ("normal".into(), "ACDEFGHIKLMNPQRSTVWY".repeat(8)),
        ]);
        let index = SeedIndex::build(&database, ChainmerParams::default()).unwrap();
        index.validate().unwrap();
        index.validate_db(&database).unwrap();
        assert_eq!(index.sequence_stats["empty-seeds"].keys, 0);
        assert!(index
            .postings
            .values()
            .flatten()
            .all(|hit| hit.id != "empty-seeds"));
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
            for seed in seeds {
                let posting = SeedPosting {
                    id: id.clone(),
                    p1: seed.p1,
                    p2: seed.p2,
                };
                assert!(index.postings[&seed.key].contains(&posting));
            }
        }
        index.validate().unwrap();
        index.validate_db(&database).unwrap();
    }

    #[test]
    fn save_overwrites_an_existing_index_and_cleans_temporary_files() {
        let path = std::env::temp_dir().join(format!(
            "chainmer-safe-save-{}-{}.json",
            std::process::id(),
            WRITE_COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
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
        assert!(index.postings.values().flatten().all(|hit| hit.p2 > hit.p1));
        assert_eq!(index.sequence_stats["a"].keys, index.postings.len());
        assert!(index.sequence_stats["a"].idf_weight.is_finite());
        assert_eq!(index.sequence_stats["a"].idf_weight, 0.0);
        let mut changed = db();
        changed.insert("b".into(), "ACDE".into());
        assert!(index.validate_db(&changed).is_err());
        let mut changed = db();
        changed.get_mut("a").unwrap().push('A');
        assert!(index.validate_db(&changed).is_err());
    }
}
