use crate::core::{orthochain_keys, OrthoParams};
use crate::fasta::read_fasta;
use crate::index::{database_keys, key_hex, parse_key_hex, SeedIndex};
use crate::search::{search, search_weighted_index};
use crate::seeds::{
    clean_protein, kmers, minimizers, spaced_kmers, strobemers, syncmers, to_class_string,
    SecPolicy,
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
    name = "orthochain",
    version,
    about = "OrthoChain v0.5.0：双通道位置种子、共线链与检索"
)]
pub struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// 统计或导出每个序列的唯一种子键
    Seed(SeedArgs),
    /// 构建带完整元数据的 OrthoChain 倒排索引
    Index(IndexArgs),
    /// 使用索引参数检索，并校验数据库内容摘要
    Search(SearchArgs),
    /// 在同一 FASTA 上比较种子数量和提取耗时
    Compare(CompareArgs),
    /// 运行明确标注的合成序列示例
    Demo,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Method {
    OrthoChain,
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
            Self::OrthoChain => "OrthoChain",
            Self::Kmer => "k-mer(k=5)",
            Self::Spaced => "spaced(11101101)",
            Self::Minimizer => "minimizer(k5w8)",
            Self::Syncmer => "syncmer(k5s2,open)",
            Self::Strobemer => "simplified-strobemer(k4)",
            Self::Ckmer => "ck-mer(k=5)",
            Self::Subseqhash2 => "SubseqHash2-protein-scalar(n25k11d5r11)",
        }
    }
}

fn extract(method: Method, sequence: &str) -> Result<BTreeSet<u128>> {
    let cleaned = || clean_protein(sequence, SecPolicy::C);
    match method {
        Method::OrthoChain => orthochain_keys(sequence, &OrthoParams::default()),
        Method::Kmer => Ok(kmers(&cleaned(), 5)?.into_iter().map(|x| x.1).collect()),
        Method::Spaced => Ok(spaced_kmers(&cleaned(), "11101101", 6)?
            .into_iter()
            .map(|x| x.1)
            .collect()),
        Method::Minimizer => minimizers(&cleaned(), 5, 8),
        Method::Syncmer => syncmers(&cleaned(), 5, 2, false),
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
    #[arg(short, long, value_enum, default_value = "orthochain")]
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
    params: OrthoCliParams,
}

#[derive(Args)]
struct SearchArgs {
    #[arg(short, long)]
    db: PathBuf,
    #[arg(short, long)]
    query: PathBuf,
    #[arg(short = 'x', long)]
    index: Option<PathBuf>,
    #[arg(long, default_value_t = 1)]
    min_shared: usize,
    #[arg(long)]
    min_score: Option<f64>,
    /// 可选 JSON：32 位十六进制键 -> 浮点权重
    #[arg(short, long)]
    weights: Option<PathBuf>,
    #[arg(long, default_value_t = 5)]
    top: usize,
}

#[derive(Args)]
struct CompareArgs {
    #[arg(short, long)]
    input: PathBuf,
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Args, Clone)]
struct OrthoCliParams {
    #[arg(long, default_value_t = 5)]
    kc: usize,
    #[arg(long, default_value_t = 8)]
    class_kc: usize,
    #[arg(long, default_value_t = 30)]
    gmin: usize,
    #[arg(long, default_value_t = 300)]
    gmax: usize,
    #[arg(long, default_value_t = 15)]
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

impl OrthoCliParams {
    fn build(self) -> OrthoParams {
        let mut p = OrthoParams {
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
            ..OrthoParams::default()
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
    let index = match args.index {
        Some(path) => {
            let idx = SeedIndex::load(path)?;
            idx.validate_db(&db)?;
            idx
        }
        None => SeedIndex::build(&db, OrthoParams::default())?,
    };
    let weights = args.weights.as_ref().map(read_weights).transpose()?;
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
        } else {
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
        }
    }
    Ok(())
}

fn cmd_compare(args: CompareArgs) -> Result<()> {
    let sequences = read_fasta(&args.input)?;
    let methods = [
        Method::OrthoChain,
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
    let index = SeedIndex::build(&db, OrthoParams::default())?;
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
    let keys = database_keys(&db, &OrthoParams::default())?;
    assert_eq!(keys.len(), 2);
    println!("NOTE: synthetic illustrative sequences, not biological evidence.");
    Ok(())
}
