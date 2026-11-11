use crate::core::ChainmerParams;
use crate::index::{database_meta, SeedIndex};
use crate::{HASH_SCHEME, INDEX_FORMAT, VERSION};
use anyhow::{Context, Result};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

/// 索引缓存选项：对应 CLI 的 `--rebuild-index` / `--cache-dir` / `--no-cache`。
#[derive(Clone, Debug, Default)]
pub struct CacheOptions {
    pub cache_dir: Option<PathBuf>,
    pub no_cache: bool,
    pub rebuild: bool,
}

/// 平台缓存目录：优先 `CHAINMER_CACHE_DIR`，其次 Windows `LOCALAPPDATA` /
/// Unix `XDG_CACHE_HOME` 或 `~/.cache`，最后退回临时目录。
fn default_cache_dir() -> PathBuf {
    if let Some(dir) = std::env::var_os("CHAINMER_CACHE_DIR") {
        return PathBuf::from(dir);
    }
    #[cfg(windows)]
    {
        if let Some(base) = std::env::var_os("LOCALAPPDATA") {
            return PathBuf::from(base).join("chainmer").join("cache");
        }
    }
    #[cfg(not(windows))]
    {
        if let Some(base) = std::env::var_os("XDG_CACHE_HOME") {
            return PathBuf::from(base).join("chainmer");
        }
        if let Some(home) = std::env::var_os("HOME") {
            return PathBuf::from(home).join(".cache").join("chainmer");
        }
    }
    std::env::temp_dir().join("chainmer-cache")
}

/// 缓存身份 = 索引格式 + 算法 + 版本 + 哈希方案 + 参数 + 数据库内容摘要与序列数。
/// 任一项变化都会得到不同的缓存文件名，因此“缓存文件存在”不等于“缓存有效”。
fn cache_filename(db: &BTreeMap<String, String>, params: &ChainmerParams) -> Result<String> {
    let db_meta = database_meta(db);
    let mut hasher = blake3::Hasher::new();
    hasher.update(INDEX_FORMAT.as_bytes());
    hasher.update(b"\0");
    hasher.update(b"Chainmer-dual-chain-xmasked");
    hasher.update(b"\0");
    hasher.update(VERSION.as_bytes());
    hasher.update(b"\0");
    hasher.update(HASH_SCHEME.as_bytes());
    hasher.update(b"\0");
    hasher.update(&(db_meta.n_sequences as u64).to_le_bytes());
    hasher.update(db_meta.content_blake3.as_bytes());
    hasher.update(b"\0");
    hasher.update(&serde_json::to_vec(params).context("无法序列化索引参数用于缓存身份")?);
    Ok(format!("chainmer-{}.ochidx", hasher.finalize().to_hex()))
}

/// 单一入口：显式 `-x` → 严格加载并校验；否则按缓存身份查找/构建。
///
/// - 显式索引与数据库不匹配时直接报错，不偷偷重建；
/// - 缓存命中且通过 `validate_db` 时直接复用，避免每次 search 重新提取全库种子；
/// - 缓存缺失/损坏/`--rebuild-index` 时构建到原子替换的缓存文件；
/// - `--no-cache` 为明确的一次性内存构建，不读也不写缓存。
pub fn load_or_build_index(
    db: &BTreeMap<String, String>,
    params: &ChainmerParams,
    explicit_index: Option<&Path>,
    options: &CacheOptions,
) -> Result<SeedIndex> {
    if let Some(path) = explicit_index {
        let index = SeedIndex::load(path)?;
        index.validate_db(db)?;
        return Ok(index);
    }

    if options.no_cache {
        let index = SeedIndex::build(db, params.clone())?;
        index.validate_db(db)?;
        eprintln!(
            "index: no-cache -> built in memory (targets={} keys={} postings={})",
            index.targets.len(),
            index.keys.len(),
            index.postings.len()
        );
        return Ok(index);
    }

    let cache_dir = options.cache_dir.clone().unwrap_or_else(default_cache_dir);
    let cache_path = cache_dir.join(cache_filename(db, params)?);

    if !options.rebuild && cache_path.is_file() {
        match SeedIndex::load(&cache_path).and_then(|index| {
            index.validate_db(db)?;
            Ok(index)
        }) {
            Ok(index) => {
                eprintln!(
                    "index: cache hit path={} (targets={} keys={} postings={})",
                    cache_path.display(),
                    index.targets.len(),
                    index.keys.len(),
                    index.postings.len()
                );
                return Ok(index);
            }
            Err(error) => {
                eprintln!(
                    "index: cache invalid ({error:#}), rebuilding path={}",
                    cache_path.display()
                );
            }
        }
    }

    let index = SeedIndex::build(db, params.clone())?;
    index.validate_db(db)?;
    fs::create_dir_all(&cache_dir)
        .with_context(|| format!("无法创建索引缓存目录：{}", cache_dir.display()))?;
    index.save(&cache_path)?;
    eprintln!(
        "index: cache miss -> built path={} (targets={} keys={} postings={})",
        cache_path.display(),
        index.targets.len(),
        index.keys.len(),
        index.postings.len()
    );
    Ok(index)
}
