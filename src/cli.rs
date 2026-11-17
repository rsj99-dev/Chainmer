use crate::adaptive::{adaptive_search, AdaptiveParams, FamilyContext};
use crate::core::{chainmer_keys, chainmer_keys_dense_full, ChainmerParams};
use crate::fasta::read_fasta;
use crate::index::{database_keys, key_hex, parse_key_hex, SeedIndex};
use crate::index_cache::{load_or_build_index, CacheOptions};
use crate::search::{search, search_adaptive, search_weighted_index};
use crate::seeds::{
    clean_protein, kmers, minimizers, spaced_kmers, strobemers, syncmers_with_mode,
    to_class_string, SecPolicy, SyncmerMode,
};
use crate::subseqhash2::{subseqhash2_keys, SubseqHash2Params};
use anyhow::{bail, Context, Result};
use clap::{Args, Parser, Subcommand, ValueEnum};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser)]
#[command(
    name = "chainmer",
    version,
    about = "Chainmer v1.0.3：双通道位置种子＋密度混合间距通道、共线链与检索；底层默认体制自适应判别（--native 回退原生）"
)]
pub struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// 统计或导出每个序列的唯一种子键
    Seed(SeedArgs),
    /// 构建带完整元数据的 Chainmer 倒排索引
    Index(IndexArgs),
    /// 使用索引参数检索，并校验数据库内容摘要；默认体制自适应判别，--native 回退原生
    Search(SearchArgs),
    /// 体制自适应判别全排序导出（暴力参照路径，FASTA 文件序平局裁决）
    Adaptive(AdaptiveArgs),
    /// 在同一 FASTA 上比较种子数量和提取耗时
    Compare(CompareArgs),
    /// 运行明确标注的合成序列示例
    Demo,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Method {
    #[value(name = "chainmer")]
    Chainmer,
    /// 密度混合候选：双通道锚点对键 ∪ 间距 k-mer 加密键（stride 2）
    #[value(name = "chainmer-dense")]
    ChainmerDense,
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
struct IndexArgs {
    #[arg(short, long)]
    db: PathBuf,
    #[arg(short, long)]
    output: PathBuf,
    #[command(flatten)]
    params: ChainmerCliParams,
}

#[derive(Args)]
struct SearchArgs {
    #[arg(short, long)]
    db: PathBuf,
    #[arg(short, long)]
    query: PathBuf,
    #[arg(short = 'x', long)]
    index: Option<PathBuf>,
    /// 强制重建持久化索引缓存
    #[arg(long)]
    rebuild_index: bool,
    /// 索引缓存目录；缺省使用平台缓存目录（可用 CHAINMER_CACHE_DIR 覆盖）
    #[arg(long)]
    cache_dir: Option<PathBuf>,
    /// 一次性运行：不读取也不写入索引缓存
    #[arg(long)]
    no_cache: bool,
    #[arg(long, default_value_t = 1)]
    min_shared: usize,
    #[arg(long)]
    min_score: Option<f64>,
    /// 可选 JSON：32 位十六进制键 -> 浮点权重
    #[arg(short, long)]
    weights: Option<PathBuf>,
    #[arg(long, default_value_t = 5)]
    top: usize,
    /// 回退原生共线链打分（关闭体制自适应判别层）
    #[arg(long)]
    native: bool,
    /// 背景扣除强度 λ（adaptive）
    #[arg(long, default_value_t = crate::adaptive::DEFAULT_LAMBDA)]
    lambda: f64,
    /// 诊断子序列通道权重 β（adaptive）
    #[arg(long, default_value_t = crate::adaptive::DEFAULT_BETA)]
    beta: f64,
    /// 家族级架构门控阈值（adaptive；最小两两余弦低于此值判分歧体制）
    #[arg(long, default_value_t = crate::adaptive::DEFAULT_DIV_THRESHOLD)]
    div_threshold: f64,
}

#[derive(Args)]
struct AdaptiveArgs {
    /// 同源集库 FASTA（含查询自身；ID 约定 species__accession）
    #[arg(short, long)]
    db: PathBuf,
    /// 查询 FASTA（ID 必须存在于库中）
    #[arg(short, long)]
    query: PathBuf,
    /// TSV 输出路径；缺省打印到 stdout
    #[arg(short, long)]
    output: Option<PathBuf>,
    /// 背景扣除强度 λ
    #[arg(long, default_value_t = crate::adaptive::DEFAULT_LAMBDA)]
    lambda: f64,
    /// 诊断子序列通道权重 β
    #[arg(long, default_value_t = crate::adaptive::DEFAULT_BETA)]
    beta: f64,
    /// 家族级架构门控阈值（最小两两余弦低于此值判分歧体制）
    #[arg(long, default_value_t = crate::adaptive::DEFAULT_DIV_THRESHOLD)]
    div_threshold: f64,
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
        Command::Index(args) => cmd_index(args),
        Command::Search(args) => cmd_search(args),
        Command::Adaptive(args) => cmd_adaptive(args),
        Command::Compare(args) => cmd_compare(args),
        Command::Demo => cmd_demo(),
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

fn cmd_index(args: IndexArgs) -> Result<()> {
    let db = read_fasta(&args.db)?;
    let index = SeedIndex::build(&db, args.params.build())?;
    let keys = index.postings.len();
    index.save(&args.output)?;
    println!(
        "wrote {}: keys={} seqs={} version={}",
        args.output.display(),
        keys,
        db.len(),
        crate::VERSION
    );
    Ok(())
}

fn read_weights(path: &PathBuf) -> Result<BTreeMap<u128, f64>> {
    let text = fs::read_to_string(path)?;
    let raw: BTreeMap<String, f64> = serde_json::from_str(&text)?;
    if raw.values().any(|w| !w.is_finite()) {
        bail!("权重必须为有限数值");
    }
    let mut weights = BTreeMap::new();
    for (key, value) in raw {
        let parsed = parse_key_hex(&key)?;
        if weights.insert(parsed, value).is_some() {
            bail!("权重文件包含重复解析键：{key}");
        }
    }
    Ok(weights)
}

fn cmd_search(args: SearchArgs) -> Result<()> {
    if args.top == 0 {
        bail!("top 必须 >= 1");
    }
    if args.min_shared == 0 {
        bail!("min_shared 必须 >= 1");
    }
    let db = read_fasta(&args.db)?;
    let index = load_or_build_index(
        &db,
        &ChainmerParams::default(),
        args.index.as_deref(),
        &CacheOptions {
            cache_dir: args.cache_dir.clone(),
            no_cache: args.no_cache,
            rebuild: args.rebuild_index,
        },
    )?;
    let weights = args.weights.as_ref().map(read_weights).transpose()?;
    let ap = AdaptiveParams {
        lambda: args.lambda,
        beta: args.beta,
        div_threshold: args.div_threshold,
    };
    if !args.native && weights.is_none() {
        eprintln!(
            "search: adaptive (default) fam_divergent={} lambda={} beta={} div_threshold={}",
            index.meta.fam_divergent, ap.lambda, ap.beta, ap.div_threshold
        );
    }
    for (qid, query) in read_fasta(&args.query)? {
        println!(">{qid}");
        if let Some(weights) = weights.as_ref() {
            let hits = search_weighted_index(
                &query,
                &index,
                &index.meta.params,
                weights,
                args.min_score,
                args.min_shared,
            )?;
            for (id, score) in hits.into_iter().take(args.top) {
                println!("{id}\t{score:.8}");
            }
        } else if args.native {
            if args.min_score.is_some() {
                bail!("min_score 仅适用于 --weights 加权检索");
            }
            let hits = search(
                &query,
                &index,
                &index.meta.params,
                args.min_shared,
                &BTreeSet::new(),
                0,
                None,
            )?;
            for hit in hits.into_iter().take(args.top) {
                println!("{}\t{:.8}", hit.id, hit.score);
            }
        } else {
            if args.min_score.is_some() {
                bail!("min_score 仅适用于 --weights 加权检索");
            }
            let result = search_adaptive(
                &query,
                &qid,
                &index,
                &index.meta.params,
                &ap,
                args.min_shared,
            )?;
            eprintln!("  regime={}", result.regime.label());
            for hit in result.hits.into_iter().take(args.top) {
                println!("{}\t{:.8}", hit.id, hit.score);
            }
        }
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

fn cmd_adaptive(args: AdaptiveArgs) -> Result<()> {
    if !args.lambda.is_finite() || !args.beta.is_finite() || !args.div_threshold.is_finite() {
        bail!("lambda、beta 与 div_threshold 必须为有限数值");
    }
    let db = read_fasta(&args.db)?;
    let ids = fasta_order(&args.db)?;
    let seqs: Vec<String> = ids
        .iter()
        .map(|id| db.get(id).cloned().context(format!("库中缺少序列：{id}")))
        .collect::<Result<_>>()?;
    let ap = AdaptiveParams {
        lambda: args.lambda,
        beta: args.beta,
        div_threshold: args.div_threshold,
    };
    let cp = ChainmerParams::default();
    let ctx = FamilyContext::build(ids, &seqs, &ap, &cp)?;
    eprintln!(
        "adaptive: seqs={} fam_divergent={} lambda={} beta={} div_threshold={}",
        ctx.ids.len(),
        ctx.fam_divergent,
        ap.lambda,
        ap.beta,
        ap.div_threshold
    );
    let mut out: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(BufWriter::new(fs::File::create(path)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };
    writeln!(out, "query\tregime\trank\tid\tscore\tgate")?;
    for (qid, qseq) in read_fasta(&args.query)? {
        if !db.contains_key(&qid) {
            bail!("查询 {qid} 不在库中；adaptive 要求查询属于同源集");
        }
        let result = adaptive_search(&ctx, &qid, &qseq, &ap, &cp)?;
        for (rank, hit) in result.hits.iter().enumerate() {
            writeln!(
                out,
                "{qid}\t{}\t{}\t{}\t{:.10}\t{:.10}",
                result.regime.label(),
                rank + 1,
                hit.id,
                hit.score,
                hit.gate
            )?;
        }
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

fn cmd_demo() -> Result<()> {
    let db=BTreeMap::from([
        ("MvhA_synthetic".into(),"MVKAVLAVLAKGGTIAAQLSPEHMPVRMRMLDNNEAAELKALGVPFGRIRLPGGQVAEPIKVDANGVITDAQSSSLIFRGFERILQHQDPRDAAFLTQRI".into()),
        ("FrhA_synthetic".into(),"MTNKKAAAWLVGAIVGGLGMAEAKLDDKAEAKKAAKEAKKAAELKALGVPFGRIRLPGGQVAEPIKVDANGVITDAQSSSLIFRGFERILQHQDPRDAA".into()),
    ]);
    let index = SeedIndex::build(&db, ChainmerParams::default())?;
    for (id, seq) in &db {
        let hits = search(
            seq,
            &index,
            &index.meta.params,
            1,
            &BTreeSet::new(),
            0,
            None,
        )?;
        println!("{id}: {hits:?}");
    }
    let keys = database_keys(&db, &ChainmerParams::default())?;
    assert_eq!(keys.len(), 2);
    println!("NOTE: synthetic illustrative sequences, not biological evidence.");
    Ok(())
}
