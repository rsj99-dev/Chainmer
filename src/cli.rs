use crate::adaptive::{family_divergent, spacing_hist, subset_keys_vec, AdaptiveParams};
use crate::core::{
    chainmer_keys, chainmer_keys_dense_full, chainmer_keys_vec, chainmer_pair_keys, ChainmerParams,
};
use crate::fasta::read_fasta;
use crate::index::key_hex;
use crate::seeds::{
    clean_protein, kmers, minimizers, spaced_kmers, strobemers, syncmers_with_mode,
    to_class_string, SecPolicy, SyncmerMode,
};
use crate::subseqhash2::{subseqhash2_keys, SubseqHash2Params};
use anyhow::{bail, Context, Result};
use clap::{Args, Parser, Subcommand, ValueEnum};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser)]
#[command(
    name = "chainmer",
    version,
    about = "Chainmer v1.1.0：双通道位置种子＋密度混合间距通道、共线链与检索；底层默认体制自适应判别（--native 回退原生）"
)]
pub struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// 统计或导出每个序列的唯一种子键
    Seed(SeedArgs),
    /// 两段式检索：k-mer 候选 + Chainmer 判别重排（实验）
    Stage2(Stage2Args),
    /// 在同一 FASTA 上比较种子数量和提取耗时
    Compare(CompareArgs),
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Method {
    #[value(name = "chainmer")]
    Chainmer,
    /// 密度混合候选：双通道锚点对键 ∪ 间距 k-mer 加密键（stride 2）
    #[value(name = "chainmer-dense")]
    ChainmerDense,
    /// 仅理化类别（CLASS）通道锚点对键（实验对照）
    #[value(name = "chainmer-class")]
    ChainmerClass,
    Kmer,
    Spaced,
    Minimizer,
    Syncmer,
    Strobemer,
    Ckmer,
    Subseqhash2,
}

impl Method {
    fn name(self) -> &'static str {
        match self {
            Self::Chainmer => "Chainmer(RAW5+CLASS8+dense)",
            Self::ChainmerDense => "Chainmer+dense(RAW5+CLASS8+spaced10-11/stride2)",
            Self::ChainmerClass => "Chainmer-class(CLASS8 anchor-pair only)",
            Self::Kmer => "k-mer(k=5)",
            Self::Spaced => "spaced(11101101)",
            Self::Minimizer => "minimizer(k5w8)",
            Self::Syncmer => "syncmer(k5s2,closed)",
            Self::Strobemer => "simplified-strobemer(k4)",
            Self::Ckmer => "ck-mer(k=5)",
            Self::Subseqhash2 => "SubseqHash2-protein-scalar(n25k11d5r11)",
        }
    }
}

fn extract(method: Method, sequence: &str) -> Result<BTreeSet<u128>> {
    let cleaned = || clean_protein(sequence, SecPolicy::C);
    match method {
        Method::Chainmer => chainmer_keys(sequence, &ChainmerParams::default()),
        Method::ChainmerClass => {
            let mut p = ChainmerParams::default();
            p.dual_track = false;
            p.use_class = true;
            p.dense_masks.clear();
            chainmer_pair_keys(sequence, &p)
        }
        Method::ChainmerDense => {
            let stride = std::env::var("CHAINMER_DENSE_STRIDE")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .unwrap_or(2);
            let class_dense = std::env::var("CHAINMER_DENSE_CLASS")
                .map(|v| v == "1")
                .unwrap_or(false);
            let mut p = ChainmerParams::default();
            if let Some(n) = std::env::var("CHAINMER_DENSE_PAIRS")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
            {
                p.max_pairs_per_anchor = n;
            }
            chainmer_keys_dense_full(sequence, &p, stride, class_dense)
        }
        Method::Kmer => Ok(kmers(&cleaned(), 5)?.into_iter().map(|x| x.1).collect()),
        Method::Spaced => Ok(spaced_kmers(&cleaned(), "11101101", 6)?
            .into_iter()
            .map(|x| x.1)
            .collect()),
        Method::Minimizer => minimizers(&cleaned(), 5, 8),
        Method::Syncmer => syncmers_with_mode(&cleaned(), 5, 2, SyncmerMode::Closed),
        Method::Strobemer => strobemers(&cleaned(), 4, 6, 30, 2),
        Method::Ckmer => Ok(kmers(&to_class_string(&cleaned()), 5)?
            .into_iter()
            .map(|x| x.1)
            .collect()),
        Method::Subseqhash2 => subseqhash2_keys(sequence, &SubseqHash2Params::default()),
    }
}

#[derive(Args)]
struct SeedArgs {
    #[arg(short, long)]
    input: PathBuf,
    #[arg(short, long, value_enum, default_value = "chainmer")]
    method: Method,
    /// 输出每个 ID 的全部键，而不是只统计数量
    #[arg(long)]
    keys: bool,
    /// TSV 输出路径；--keys 时必填
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Args)]
struct CompareArgs {
    #[arg(short, long)]
    input: PathBuf,
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Args, Clone)]
struct ChainmerCliParams {
    #[arg(long, default_value_t = 5)]
    kc: usize,
    #[arg(long, default_value_t = 8)]
    class_kc: usize,
    #[arg(long, default_value_t = 30)]
    gmin: usize,
    #[arg(long, default_value_t = 300)]
    gmax: usize,
    #[arg(long, default_value_t = 5)]
    delta: usize,
    #[arg(long, default_value_t = 6)]
    w1: usize,
    #[arg(long, default_value_t = 12)]
    max_pairs_per_anchor: usize,
    #[arg(long, default_value_t = 30)]
    chain_band: usize,
    /// 只使用原始氨基酸通道
    #[arg(long, conflicts_with = "class_only")]
    raw: bool,
    /// 只使用替换等价类通道，用于消融
    #[arg(long)]
    class_only: bool,
    /// 同时关闭距离过滤和间距分桶
    #[arg(long)]
    no_gap: bool,
    #[arg(long)]
    no_gap_filter: bool,
    #[arg(long)]
    no_gap_bucket: bool,
    #[arg(long, default_value = "C")]
    sec_policy: String,
}

impl ChainmerCliParams {
    fn build(self) -> ChainmerParams {
        let mut p = ChainmerParams {
            kc: self.kc,
            class_kc: self.class_kc,
            gmin: self.gmin,
            gmax: self.gmax,
            delta: self.delta,
            w1: self.w1,
            max_pairs_per_anchor: self.max_pairs_per_anchor,
            chain_band: self.chain_band,
            use_class: !self.raw,
            dual_track: !self.raw && !self.class_only,
            sec_policy: self.sec_policy,
            ..ChainmerParams::default()
        };
        if self.no_gap {
            p.set_use_gap(false);
        }
        if self.no_gap_filter {
            p.gap_filter = false;
        }
        if self.no_gap_bucket {
            p.gap_bucket = false;
        }
        p
    }
}

pub fn run() -> Result<()> {
    match Cli::parse().command {
        Command::Seed(args) => cmd_seed(args),
        Command::Stage2(args) => cmd_stage2(args),
        Command::Compare(args) => cmd_compare(args),
    }
}

fn cmd_seed(args: SeedArgs) -> Result<()> {
    let sequences = read_fasta(&args.input)?;
    if args.keys {
        let path = args
            .output
            .as_ref()
            .context("--keys 需要同时指定 --output")?;
        let mut writer = BufWriter::new(fs::File::create(path)?);
        writeln!(writer, "# method={}", args.method.name())?;
        for (id, seq) in &sequences {
            for key in extract(args.method, seq)? {
                writeln!(writer, "{id}\t{}", key_hex(key))?;
            }
        }
        writer.flush()?;
        eprintln!(
            "已写入 {}：{} 条序列，方法 {}",
            path.display(),
            sequences.len(),
            args.method.name()
        );
    } else {
        let mut total = 0usize;
        for (id, seq) in &sequences {
            let count = extract(args.method, seq)?.len();
            total += count;
            println!("{id}\t{count}");
        }
        println!(
            "# sequences={} total_keys={} method={}",
            sequences.len(),
            total,
            args.method.name()
        );
    }
    Ok(())
}

/// 按文件顺序提取 FASTA ID（平局裁决与物种分组顺序依赖文件序）。
fn fasta_order(path: &PathBuf) -> Result<Vec<String>> {
    let text = fs::read_to_string(path)?;
    Ok(text
        .lines()
        .filter_map(|line| line.strip_prefix('>'))
        .filter_map(|header| header.split_whitespace().next())
        .map(|id| id.to_string())
        .collect())
}

#[derive(Args)]
struct Stage2Args {
    #[arg(short, long)]
    db: PathBuf,
    #[arg(short, long)]
    query: PathBuf,
    /// 第一段 k-mer 召回候选数
    #[arg(long, default_value_t = 200)]
    n: usize,
    /// 线程数（预计算与查询）
    #[arg(long, default_value_t = 4)]
    threads: usize,
    #[arg(long, default_value_t = 20)]
    top: usize,
    #[arg(long, default_value_t = 1.0)]
    lambda: f64,
    #[arg(long, default_value_t = 1.0)]
    beta: f64,
    #[arg(long, default_value_t = 0.8)]
    div_threshold: f64,
    #[arg(short, long)]
    output: Option<PathBuf>,
}

fn df_counts(keys: &[BTreeSet<u128>]) -> HashMap<u128, u32> {
    let mut df: HashMap<u128, u32> = HashMap::new();
    for ks in keys {
        for k in ks {
            *df.entry(*k).or_insert(0) += 1;
        }
    }
    df
}

/// IDF 加权 Jaccard（df 已预计算，O(|q∪t|)）。
fn wj(q: &BTreeSet<u128>, t: &BTreeSet<u128>, df: &HashMap<u128, u32>, n: usize) -> f64 {
    if q.is_empty() || t.is_empty() {
        return 0.0;
    }
    let w = |k: u128| (1.0 + n as f64 / (*df.get(&k).unwrap_or(&1) as f64)).ln();
    let mut inter = 0.0;
    let mut uni = 0.0;
    let (mut iq, mut it) = (q.iter().peekable(), t.iter().peekable());
    loop {
        match (iq.peek(), it.peek()) {
            (Some(a), Some(b)) => {
                if a < b {
                    uni += w(**a);
                    iq.next();
                } else if a > b {
                    uni += w(**b);
                    it.next();
                } else {
                    let ww = w(**a);
                    inter += ww;
                    uni += ww;
                    iq.next();
                    it.next();
                }
            }
            (Some(a), None) => {
                uni += w(**a);
                iq.next();
            }
            (None, Some(b)) => {
                uni += w(**b);
                it.next();
            }
            (None, None) => break,
        }
    }
    if uni > 0.0 {
        inter / uni
    } else {
        0.0
    }
}

fn w_of(k: u128, df: &HashMap<u128, u32>, n: usize) -> f64 {
    (1.0 + n as f64 / (*df.get(&k).unwrap_or(&1) as f64)).ln()
}

fn df_counts_vec(keys: &[Vec<u128>]) -> HashMap<u128, u32> {
    let mut df: HashMap<u128, u32> = HashMap::new();
    for ks in keys {
        for k in ks {
            *df.entry(*k).or_insert(0) += 1;
        }
    }
    df
}

/// 紧凑权重表：升序键数组（u64 截断键）+ f32 权重，二分查找（省内存）。
struct WTab {
    keys: Vec<u64>,
    vals: Vec<f32>,
}

impl WTab {
    fn new(df: &HashMap<u64, u32>, n: usize) -> Self {
        let mut v: Vec<(u64, f32)> = df
            .iter()
            .map(|(&k, &c)| (k, (1.0 + n as f64 / c as f64).ln() as f32))
            .collect();
        v.sort_unstable_by_key(|&(k, _)| k);
        Self {
            keys: v.iter().map(|x| x.0).collect(),
            vals: v.iter().map(|x| x.1).collect(),
        }
    }

    fn get(&self, k: &u64) -> f64 {
        match self.keys.binary_search(k) {
            Ok(i) => self.vals[i] as f64,
            Err(_) => 0.0,
        }
    }

    /// 从已排序的 (key, seq) 倒排一次扫描生成（免 df HashMap 与二次排序）。
    fn from_sorted_inv(inv: &[(u64, u32)], n: usize) -> Self {
        let mut keys: Vec<u64> = Vec::new();
        let mut vals: Vec<f32> = Vec::new();
        let mut i = 0usize;
        while i < inv.len() {
            let k = inv[i].0;
            let mut j = i + 1;
            while j < inv.len() && inv[j].0 == k {
                j += 1;
            }
            let df = (j - i) as f64;
            keys.push(k);
            vals.push((1.0 + n as f64 / df).ln() as f32);
            i = j;
        }
        Self { keys, vals }
    }
}

fn df_counts_u64(keys: &[Vec<u64>]) -> HashMap<u64, u32> {
    let mut df: HashMap<u64, u32> = HashMap::new();
    for ks in keys {
        for k in ks {
            *df.entry(*k).or_insert(0) += 1;
        }
    }
    df
}

/// u128 键截断为 u64（碰撞期望 ~1e-7）并重排去重。
fn u64_keys(v: &[u128]) -> Vec<u64> {
    let mut out: Vec<u64> = v.iter().map(|k| *k as u64).collect();
    out.sort_unstable();
    out.dedup();
    out
}

fn weight_table(df: &HashMap<u128, u32>, n: usize) -> HashMap<u128, f64> {
    df.iter()
        .map(|(&k, &c)| (k, (1.0 + n as f64 / c as f64).ln()))
        .collect()
}

/// 排序 Vec 双指针加权 Jaccard（权重查表，无分配）。
fn wj_sorted(q: &[u128], t: &[u128], w: &HashMap<u128, f64>) -> f64 {
    if q.is_empty() || t.is_empty() {
        return 0.0;
    }
    let (mut i, mut j) = (0usize, 0usize);
    let (mut inter, mut uni) = (0.0f64, 0.0f64);
    while i < q.len() && j < t.len() {
        if q[i] < t[j] {
            uni += w.get(&q[i]).copied().unwrap_or(0.0);
            i += 1;
        } else if q[i] > t[j] {
            uni += w.get(&t[j]).copied().unwrap_or(0.0);
            j += 1;
        } else {
            let ww = w.get(&q[i]).copied().unwrap_or(0.0);
            inter += ww;
            uni += ww;
            i += 1;
            j += 1;
        }
    }
    while i < q.len() {
        uni += w.get(&q[i]).copied().unwrap_or(0.0);
        i += 1;
    }
    while j < t.len() {
        uni += w.get(&t[j]).copied().unwrap_or(0.0);
        j += 1;
    }
    if uni > 0.0 {
        inter / uni
    } else {
        0.0
    }
}

fn merge_sorted(sets: &[&Vec<u64>]) -> Vec<u64> {
    let mut all: Vec<u64> = sets.iter().flat_map(|s| s.iter().copied()).collect();
    all.sort_unstable();
    all.dedup();
    all
}

/// 倒排稀疏交集：只遍历查询键的 posting（二分定位键块）。
/// 返回（每个库序列的加权交集，查询键权重和）。
fn sparse_overlap(q: &[u64], inv: &[(u64, u32)], w: &WTab, n: usize) -> (Vec<f64>, f64) {
    let mut inter = vec![0.0f64; n];
    let mut sum_q = 0.0f64;
    for k in q {
        let ww = w.get(k);
        sum_q += ww;
        let start = inv.partition_point(|&(key, _)| key < *k);
        let mut i = start;
        while i < inv.len() && inv[i].0 == *k {
            inter[inv[i].1 as usize] += ww;
            i += 1;
        }
    }
    (inter, sum_q)
}

fn par_map<T, F>(seqs: &[String], threads: usize, f: F) -> Result<Vec<T>>
where
    T: Send,
    F: Fn(&String) -> Result<T> + Sync,
{
    let nthreads = threads.max(1).min(seqs.len().max(1));
    let per = ((seqs.len() + nthreads - 1) / nthreads).max(1);
    let mut parts: Vec<(usize, Vec<T>)> = Vec::new();
    std::thread::scope(|s| -> Result<()> {
        let mut hs = Vec::new();
        for (ci, chunk) in seqs.chunks(per).enumerate() {
            let fref = &f;
            let base = ci * per;
            hs.push(s.spawn(move || -> Result<(usize, Vec<T>)> {
                let mut v = Vec::with_capacity(chunk.len());
                for x in chunk {
                    v.push(fref(x)?);
                }
                Ok((base, v))
            }));
        }
        for h in hs {
            parts.push(h.join().unwrap()?);
        }
        Ok(())
    })?;
    parts.sort_by_key(|x| x.0);
    let mut out = Vec::new();
    for (_, v) in parts {
        out.extend(v);
    }
    Ok(out)
}

fn par_ck_sk(
    seqs: &[String],
    cp: &ChainmerParams,
    threads: usize,
) -> Result<(Vec<Vec<u128>>, Vec<Vec<u128>>)> {
    let t0 = std::time::Instant::now();
    let ck = par_map(seqs, threads, |x| chainmer_keys_vec(x, cp))?;
    eprintln!("stage2: ck {:.2}s", t0.elapsed().as_secs_f64());
    let t1 = std::time::Instant::now();
    let sk = par_map(seqs, threads, |x| Ok(subset_keys_vec(x)))?;
    eprintln!("stage2: sk {:.2}s", t1.elapsed().as_secs_f64());
    Ok((ck, sk))
}

fn stage2_cache_path(db: &std::path::Path) -> PathBuf {
    let name = format!(
        "{}.stage2cache",
        db.file_name().and_then(|s| s.to_str()).unwrap_or("db")
    );
    db.with_file_name(name)
}

fn build_inv(keys: &[Vec<u64>]) -> Vec<(u64, u32)> {
    let mut inv: Vec<(u64, u32)> = Vec::new();
    for (i, ks) in keys.iter().enumerate() {
        for k in ks {
            inv.push((*k, i as u32));
        }
    }
    inv.sort_unstable_by_key(|&(k, _)| k);
    inv
}

fn save_stage2_cache(
    path: &std::path::Path,
    db_len: u64,
    db_mtime: u64,
    ck: &[Vec<u64>],
    sk: &[Vec<u64>],
    inv_ck: &[(u64, u32)],
    inv_sk: &[(u64, u32)],
) -> Result<()> {
    let mut f = BufWriter::new(fs::File::create(path)?);
    f.write_all(b"S2C2")?;
    f.write_all(&db_len.to_le_bytes())?;
    f.write_all(&db_mtime.to_le_bytes())?;
    f.write_all(&(ck.len() as u32).to_le_bytes())?;
    for arr in [ck, sk] {
        for ks in arr {
            f.write_all(&(ks.len() as u32).to_le_bytes())?;
            for k in ks {
                f.write_all(&k.to_le_bytes())?;
            }
        }
    }
    for inv in [inv_ck, inv_sk] {
        f.write_all(&(inv.len() as u64).to_le_bytes())?;
        for &(k, s) in inv {
            f.write_all(&k.to_le_bytes())?;
            f.write_all(&s.to_le_bytes())?;
        }
    }
    f.flush()?;
    Ok(())
}

fn load_stage2_cache(
    path: &std::path::Path,
    db_len: u64,
    db_mtime: u64,
    n: usize,
) -> Option<(Vec<Vec<u64>>, Vec<Vec<u64>>, Vec<(u64, u32)>, Vec<(u64, u32)>)> {
    let mut f = BufReader::new(fs::File::open(path).ok()?);
    let mut magic = [0u8; 4];
    f.read_exact(&mut magic).ok()?;
    if &magic != b"S2C2" {
        return None;
    }
    let mut b8 = [0u8; 8];
    let mut b4 = [0u8; 4];
    f.read_exact(&mut b8).ok()?;
    if u64::from_le_bytes(b8) != db_len {
        return None;
    }
    f.read_exact(&mut b8).ok()?;
    if u64::from_le_bytes(b8).abs_diff(db_mtime) > 2 {
        return None;
    }
    f.read_exact(&mut b4).ok()?;
    if u32::from_le_bytes(b4) as usize != n {
        return None;
    }
    let mut read_keys = |f: &mut BufReader<fs::File>| -> Option<Vec<Vec<u64>>> {
        let mut out: Vec<Vec<u64>> = Vec::with_capacity(n);
        for _ in 0..n {
            let mut b4 = [0u8; 4];
            f.read_exact(&mut b4).ok()?;
            let len = u32::from_le_bytes(b4) as usize;
            let mut v: Vec<u64> = Vec::with_capacity(len);
            for _ in 0..len {
                let mut b8 = [0u8; 8];
                f.read_exact(&mut b8).ok()?;
                v.push(u64::from_le_bytes(b8));
            }
            out.push(v);
        }
        Some(out)
    };
    let ck = read_keys(&mut f)?;
    let sk = read_keys(&mut f)?;
    let mut read_inv = |f: &mut BufReader<fs::File>| -> Option<Vec<(u64, u32)>> {
        let mut b8 = [0u8; 8];
        f.read_exact(&mut b8).ok()?;
        let len = u64::from_le_bytes(b8) as usize;
        let mut v: Vec<(u64, u32)> = Vec::with_capacity(len);
        for _ in 0..len {
            let mut b8 = [0u8; 8];
            f.read_exact(&mut b8).ok()?;
            let k = u64::from_le_bytes(b8);
            let mut b4 = [0u8; 4];
            f.read_exact(&mut b4).ok()?;
            v.push((k, u32::from_le_bytes(b4)));
        }
        Some(v)
    };
    let inv_ck = read_inv(&mut f)?;
    let inv_sk = read_inv(&mut f)?;
    Some((ck, sk, inv_ck, inv_sk))
}

fn cmd_stage2(args: Stage2Args) -> Result<()> {
    let db = read_fasta(&args.db)?;
    let ids = fasta_order(&args.db)?;
    let seqs: Vec<String> = ids
        .iter()
        .map(|id| db.get(id).cloned().context(format!("库中缺少序列：{id}")))
        .collect::<Result<_>>()?;
    let cp = ChainmerParams::default();
    let ap = AdaptiveParams {
        lambda: args.lambda,
        beta: args.beta,
        div_threshold: args.div_threshold,
    };
    let n = seqs.len();
    let n_cand = args.n.max(1);
    let t_build = std::time::Instant::now();
    let hists: Vec<Vec<f64>> = seqs.iter().map(|s| spacing_hist(s)).collect();
    let fam_divergent = family_divergent(&hists, ap.div_threshold);
    // 库级键表（升序 Vec，双指针用）与权重表（一次）
    let k3: Vec<Vec<u64>> = seqs
        .iter()
        .map(|s| {
            let v: Vec<u128> = kmers(&clean_protein(s, SecPolicy::C), 3)?
                .into_iter()
                .map(|x| x.1)
                .collect();
            Ok(u64_keys(&v))
        })
        .collect::<Result<_>>()?;
    let t_cache = std::time::Instant::now();
    let meta = fs::metadata(&args.db)?;
    let db_len = meta.len();
    let db_mtime = meta
        .modified()
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let cache_path = stage2_cache_path(&args.db);
    let (ck, sk, inv_ck, inv_sk) = match load_stage2_cache(&cache_path, db_len, db_mtime, n) {
        Some(v) => {
            eprintln!("stage2: cache hit {:.2}s", t_cache.elapsed().as_secs_f64());
            v
        }
        None => {
            let (ck128, sk128) = par_ck_sk(&seqs, &cp, args.threads)?;
            let ck: Vec<Vec<u64>> = ck128.iter().map(|v| u64_keys(v)).collect();
            let sk: Vec<Vec<u64>> = sk128.iter().map(|v| u64_keys(v)).collect();
            let inv_ck = build_inv(&ck);
            let inv_sk = build_inv(&sk);
            if let Err(e) =
                save_stage2_cache(&cache_path, db_len, db_mtime, &ck, &sk, &inv_ck, &inv_sk)
            {
                eprintln!("stage2: cache save failed: {e}");
            }
            eprintln!("stage2: built+cached {:.2}s", t_cache.elapsed().as_secs_f64());
            (ck, sk, inv_ck, inv_sk)
        }
    };
    let w_k3 = WTab::new(&df_counts_u64(&k3), n);
    let t_w = std::time::Instant::now();
    let w_ck = WTab::from_sorted_inv(&inv_ck, n);
    let w_sk = WTab::from_sorted_inv(&inv_sk, n);
    eprintln!("stage2: wtab {:.2}s", t_w.elapsed().as_secs_f64());
    let mut sum_ck: Vec<f64> = vec![0.0; n];
    {
        let mut wi = 0usize;
        for &(k, s) in &inv_ck {
            while wi < w_ck.keys.len() && w_ck.keys[wi] < k {
                wi += 1;
            }
            if wi < w_ck.keys.len() && w_ck.keys[wi] == k {
                sum_ck[s as usize] += w_ck.vals[wi] as f64;
            }
        }
    }
    let mut sum_sk: Vec<f64> = vec![0.0; n];
    {
        let mut wi = 0usize;
        for &(k, s) in &inv_sk {
            while wi < w_sk.keys.len() && w_sk.keys[wi] < k {
                wi += 1;
            }
            if wi < w_sk.keys.len() && w_sk.keys[wi] == k {
                sum_sk[s as usize] += w_sk.vals[wi] as f64;
            }
        }
    }
    // k=3 倒排（键 -> 序列索引）与每序列权重和
    let mut inv_k3: HashMap<u64, Vec<u32>> = HashMap::new();
    for (i, ks) in k3.iter().enumerate() {
        for k in ks {
            inv_k3.entry(*k).or_default().push(i as u32);
        }
    }
    let sum_k3: Vec<f64> = k3
        .iter()
        .map(|t| t.iter().map(|k| w_k3.get(k)).sum())
        .collect();
    let sp_of: Vec<String> = ids
        .iter()
        .map(|id| crate::adaptive::species_of(id).to_string())
        .collect();
    eprintln!("stage2: 预计算 {:.2}s（库 {} 条）", t_build.elapsed().as_secs_f64(), n);
    let mut out: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(BufWriter::new(fs::File::create(path)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };
    writeln!(out, "query\trank\tid\tscore\tgate")?;
    let queries: Vec<(String, String)> = read_fasta(&args.query)?.into_iter().collect();
    let process = |qid: &str, qseq: &str| -> Result<String> {
        use std::fmt::Write as _;
        let mut buf = String::new();
        let qk3v: Vec<u128> = kmers(&clean_protein(qseq, SecPolicy::C), 3)?
            .into_iter()
            .map(|x| x.1)
            .collect();
        let qk3: Vec<u64> = u64_keys(&qk3v);
        let qk: Vec<u64> = u64_keys(&chainmer_keys_vec(qseq, &cp)?);
        let qsk: Vec<u64> = u64_keys(&subset_keys_vec(qseq));
        // 第一段：k=3 倒排（稀疏命中）
        let mut inter = vec![0.0f64; n];
        let sum_q: f64 = qk3.iter().map(|k| w_k3.get(k)).sum();
        for k in &qk3 {
            if let Some(list) = inv_k3.get(k) {
                let wv = w_k3.get(k);
                for &i in list {
                    inter[i as usize] += wv;
                }
            }
        }
        let s_k: Vec<f64> = (0..n)
            .map(|i| {
                let uni = sum_q + sum_k3[i] - inter[i];
                if uni > 0.0 {
                    inter[i] / uni
                } else {
                    0.0
                }
            })
            .collect();
        let mut order: Vec<usize> = (0..n).filter(|&i| ids[i] != qid).collect();
        order.sort_by(|&a, &b| s_k[b].total_cmp(&s_k[a]));
        let cand: Vec<usize> = order
            .into_iter()
            .filter(|&i| s_k[i] > 0.0)
            .take(n_cand)
            .collect();
        // 第二段：倒排稀疏交集（只碰真正命中的序列）
        let (inter_c, sum_qc) = sparse_overlap(&qk, &inv_ck, &w_ck, n);
        let (inter_s, sum_qs) = sparse_overlap(&qsk, &inv_sk, &w_sk, n);
        let jc: Vec<f64> = cand
            .iter()
            .map(|&i| {
                let uni = sum_qc + sum_ck[i] - inter_c[i];
                if uni > 0.0 {
                    inter_c[i] / uni
                } else {
                    0.0
                }
            })
            .collect();
        let js: Vec<f64> = cand
            .iter()
            .map(|&i| {
                let uni = sum_qs + sum_sk[i] - inter_s[i];
                if uni > 0.0 {
                    inter_s[i] / uni
                } else {
                    0.0
                }
            })
            .collect();
        // 背景：同物种其余拷贝（合并有序 Vec）
        let qsp = crate::adaptive::species_of(qid);
        let mut bg_sets: Vec<&Vec<u64>> = Vec::new();
        let mut par_j = 0.0f64;
        for (i, sp) in sp_of.iter().enumerate() {
            if ids[i] == qid || sp != qsp {
                continue;
            }
            bg_sets.push(&ck[i]);
            if let Some(pos) = cand.iter().position(|&c| c == i) {
                par_j = par_j.max(jc[pos]);
            }
        }
        let bg = merge_sorted(&bg_sets);
        let jbg: Vec<f64> = if bg.is_empty() {
            vec![0.0; cand.len()]
        } else {
            let (inter_b, sum_qb) = sparse_overlap(&bg, &inv_ck, &w_ck, n);
            cand.iter()
                .map(|&i| {
                    let uni = sum_qb + sum_ck[i] - inter_b[i];
                    if uni > 0.0 {
                        inter_b[i] / uni
                    } else {
                        0.0
                    }
                })
                .collect()
        };
        let top1_j = jc.iter().cloned().fold(0.0f64, f64::max);
        let conserved = top1_j > 0.0 && par_j / top1_j > 0.7;
        let scores: Vec<f64> = (0..cand.len())
            .map(|k| {
                let base = if fam_divergent {
                    jc[k]
                } else if conserved {
                    jc[k] + ap.beta * js[k] - ap.lambda * jbg[k]
                } else {
                    jc[k] + ap.beta * js[k]
                };
                base + 0.02 * s_k[cand[k]]
            })
            .collect();
        let mut oc: Vec<usize> = (0..cand.len()).collect();
        oc.sort_by(|&a, &b| scores[b].total_cmp(&scores[a]));
        for (rank, &k) in oc.iter().take(args.top).enumerate() {
            writeln!(&mut buf, "{qid}\t{}\t{}\t{:.10}\t{:.10}", rank + 1, ids[cand[k]], scores[k], jc[k])?;
        }
        Ok(buf)
    };
    let nthreads = args.threads.max(1).min(queries.len().max(1));
    let per = ((queries.len() + nthreads - 1) / nthreads).max(1);
    let mut parts: Vec<(usize, String)> = Vec::new();
    std::thread::scope(|s| -> Result<()> {
        let mut hs = Vec::new();
        for (ci, chunk) in queries.chunks(per).enumerate() {
            let p = &process;
            let base = ci * per;
            hs.push(s.spawn(move || -> Result<(usize, String)> {
                let mut buf = String::new();
                for (qid, qseq) in chunk {
                    buf.push_str(&p(qid, qseq)?);
                }
                Ok((base, buf))
            }));
        }
        for h in hs {
            parts.push(h.join().unwrap()?);
        }
        Ok(())
    })?;
    parts.sort_by_key(|x| x.0);
    for (_, txt) in parts {
        out.write_all(txt.as_bytes())?;
    }
    out.flush()?;
    if let Some(path) = &args.output {
        eprintln!("已写入 {}", path.display());
    }
    Ok(())
}

fn cmd_compare(args: CompareArgs) -> Result<()> {
    let sequences = read_fasta(&args.input)?;
    let methods = [
        Method::Chainmer,
        Method::Subseqhash2,
        Method::Kmer,
        Method::Spaced,
        Method::Minimizer,
        Method::Syncmer,
        Method::Strobemer,
        Method::Ckmer,
    ];
    let mut text =
        String::from("method\tsequences\ttotal_keys\tmean_keys_per_seq\textract_time_ms\n");
    for method in methods {
        let start = Instant::now();
        let mut total = 0usize;
        for sequence in sequences.values() {
            total += extract(method, sequence)?.len();
        }
        let elapsed = start.elapsed().as_secs_f64() * 1000.0;
        text.push_str(&format!(
            "{}\t{}\t{}\t{:.3}\t{:.3}\n",
            method.name(),
            sequences.len(),
            total,
            total as f64 / sequences.len() as f64,
            elapsed
        ));
    }
    if let Some(path) = args.output {
        fs::write(&path, &text)?;
        println!("wrote {}", path.display());
    } else {
        print!("{text}");
    }
    Ok(())
}
