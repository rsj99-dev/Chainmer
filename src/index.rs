use crate::core::{evoseed_seeds, EvoParams, PositionedSeed};
use crate::seeds::class_map;
use crate::{HASH_SCHEME, INDEX_FORMAT, VERSION};
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

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
    pub params: EvoParams,
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
    params: &EvoParams,
) -> Result<BTreeMap<String, BTreeSet<PositionedSeed>>> {
    sequences
        .iter()
        .map(|(id, seq)| Ok((id.clone(), evoseed_seeds(seq, params)?)))
        .collect()
}

pub fn database_keys(
    sequences: &BTreeMap<String, String>,
    params: &EvoParams,
) -> Result<BTreeMap<String, BTreeSet<u128>>> {
    Ok(database_seeds(sequences, params)?
        .into_iter()
        .map(|(id, seeds)| (id, seeds.into_iter().map(|seed| seed.key).collect()))
        .collect())
}

impl SeedIndex {
    pub fn build(sequences: &BTreeMap<String, String>, params: EvoParams) -> Result<Self> {
        params.validate()?;
        if sequences.is_empty() {
            bail!("空数据库不能建立索引");
        }
        let database_seeds = database_seeds(sequences, &params)?;
        let mut postings: BTreeMap<u128, Vec<SeedPosting>> = BTreeMap::new();
        let mut target_keys: BTreeMap<String, BTreeSet<u128>> = BTreeMap::new();
        for (id, seeds) in database_seeds {
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
            ((n + 1) as f64 / (df + 1) as f64).ln()
        };
        let sequence_stats = target_keys
            .into_iter()
            .map(|(id, keys)| {
                let stat = SequenceStats {
                    keys: keys.len(),
                    idf_weight: keys.iter().map(|key| idf(key)).sum(),
                };
                (id, stat)
            })
            .collect();
        Ok(Self {
            meta: IndexMeta {
                format: INDEX_FORMAT.into(),
                algorithm: "EvoSeed-dual-chain".into(),
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
        ((self.meta.db.n_sequences + 1) as f64 / (df + 1) as f64).ln()
    }

    pub fn validate_db(&self, sequences: &BTreeMap<String, String>) -> Result<()> {
        let actual = database_meta(sequences);
        if self.meta.db != actual {
            bail!("数据库与索引不匹配：序列 ID、顺序无关内容摘要或序列数不同");
        }
        Ok(())
    }

    pub fn validate(&self) -> Result<()> {
        if self.meta.format != INDEX_FORMAT {
            bail!("不支持的索引格式：{}", self.meta.format);
        }
        if self.meta.algorithm != "EvoSeed-dual-chain" || self.meta.version != VERSION {
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
        if self.sequence_stats.len() != self.meta.db.n_sequences {
            bail!("索引序列统计数量与数据库元数据不一致");
        }
        if self
            .sequence_stats
            .values()
            .any(|stat| !stat.idf_weight.is_finite() || stat.idf_weight < 0.0)
        {
            bail!("索引包含无效的序列 IDF 权重");
        }
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
            if self.document_freq.get(key) != Some(&df) {
                bail!("索引键频率与 postings 不一致");
            }
        }
        if self.document_freq.len() != self.postings.len() {
            bail!("索引键频率表与 postings 大小不一致");
        }
        Ok(())
    }

    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        self.validate()?;
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
        fs::write(path.as_ref(), text)
            .with_context(|| format!("无法写入索引：{}", path.as_ref().display()))?;
        Ok(())
    }

    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let text = fs::read_to_string(path.as_ref())?;
        let file: IndexFile =
            serde_json::from_str(&text).context("索引必须符合 evoseed-index/3；旧版索引需重建")?;
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
    fn metadata_database_digest_and_positional_stats_are_real() {
        let mut params = EvoParams::default();
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
