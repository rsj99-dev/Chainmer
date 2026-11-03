use crate::core::{evoseed_keys, EvoParams};
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

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeedIndex {
    pub meta: IndexMeta,
    pub postings: BTreeMap<u128, Vec<String>>,
}

#[derive(Serialize, Deserialize)]
struct IndexFile {
    meta: IndexMeta,
    index: BTreeMap<String, Vec<String>>,
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

pub fn database_keys(
    sequences: &BTreeMap<String, String>,
    params: &EvoParams,
) -> Result<BTreeMap<String, BTreeSet<u128>>> {
    sequences
        .iter()
        .map(|(id, seq)| Ok((id.clone(), evoseed_keys(seq, params)?)))
        .collect()
}

impl SeedIndex {
    pub fn build(sequences: &BTreeMap<String, String>, params: EvoParams) -> Result<Self> {
        params.validate()?;
        if sequences.is_empty() {
            bail!("空数据库不能建立索引");
        }
        let mut postings: BTreeMap<u128, Vec<String>> = BTreeMap::new();
        for (id, keys) in database_keys(sequences, &params)? {
            for key in keys {
                postings.entry(key).or_default().push(id.clone());
            }
        }
        Ok(Self {
            meta: IndexMeta {
                format: INDEX_FORMAT.into(),
                algorithm: "EvoSeed".into(),
                version: VERSION.into(),
                params,
                class_map: class_map(),
                hash_scheme: HASH_SCHEME.into(),
                db: database_meta(sequences),
            },
            postings,
        })
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
        if self.meta.algorithm != "EvoSeed" || self.meta.version != VERSION {
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
        for ids in self.postings.values() {
            if ids.is_empty() {
                bail!("索引包含空 postings");
            }
            let unique: BTreeSet<_> = ids.iter().collect();
            if unique.len() != ids.len() {
                bail!("索引 postings 包含重复序列 ID");
            }
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
                .map(|(key, ids)| {
                    let mut ids = ids.clone();
                    ids.sort();
                    (key_hex(*key), ids)
                })
                .collect(),
        };
        let text = serde_json::to_string_pretty(&file)?;
        fs::write(path.as_ref(), text)
            .with_context(|| format!("无法写入索引：{}", path.as_ref().display()))?;
        Ok(())
    }

    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let text = fs::read_to_string(path.as_ref())
            .with_context(|| format!("无法读取索引：{}", path.as_ref().display()))?;
        let file: IndexFile = serde_json::from_str(&text)
            .context("索引必须符合 evoseed-index/2；旧版无元数据索引需重建")?;
        let mut postings = BTreeMap::new();
        for (key, ids) in file.index {
            let parsed = parse_key_hex(&key)?;
            if postings.insert(parsed, ids).is_some() {
                bail!("索引包含重复解析键：{key}");
            }
        }
        let index = Self {
            meta: file.meta,
            postings,
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
    fn metadata_and_database_digest_are_real() {
        let mut params = EvoParams::default();
        params.set_use_gap(false);
        let index = SeedIndex::build(&db(), params).unwrap();
        assert!(!index.meta.params.gap_filter && !index.meta.params.gap_bucket);
        assert!(!index.meta.class_map.is_empty());
        let mut changed = db();
        changed.insert("b".into(), "ACDE".into());
        assert!(index.validate_db(&changed).is_err());
        let mut changed = db();
        changed.get_mut("a").unwrap().push('A');
        assert!(index.validate_db(&changed).is_err());
    }
}
