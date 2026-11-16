//! seed-baselines 导出 CLI：与论文 tools/seed_dump.rs 的六类基线分支逐字节兼容。
//!
//! 输出格式：sequence_id <TAB> key_hex <TAB> pos
//! pos 为该键首次出现的位置（集合型种子记 -1）。统计信息写入 stderr。

use anyhow::{bail, Result};
use seed_baselines::fasta::read_fasta;
use seed_baselines::seeds::{
    clean_protein, kmers, minimizers, spaced_kmers, strobemers, syncmers_with_mode, SecPolicy,
    SyncmerMode,
};
use seed_baselines::subseqhash2::{subseqhash2_keys, SubseqHash2Params};
use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

struct Args {
    input: String,
    output: String,
    method: String,
    k: usize,
    w: usize,
    s: usize,
    w1: usize,
    w2: usize,
    t: usize,
    stride: usize,
    window: usize,
    subseq: usize,
    d: usize,
    repeats: usize,
    masks: String,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: String::new(),
            output: String::new(),
            method: String::from("kmer"),
            k: 7,
            w: 1,
            s: 3,
            w1: 6,
            w2: 30,
            t: 2,
            stride: 1,
            window: 25,
            subseq: 11,
            d: 5,
            repeats: 5,
            masks: String::from("11101101"),
        }
    }
}

fn parse_args() -> Result<Args> {
    let mut args = Args::default();
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        let key = argv[i].trim_start_matches("--").to_string();
        let value = argv
            .get(i + 1)
            .cloned()
            .ok_or_else(|| anyhow::anyhow!("参数 {key} 缺少取值"))?;
        match key.as_str() {
            "input" => args.input = value,
            "output" => args.output = value,
            "method" => args.method = value,
            "k" => args.k = value.parse()?,
            "w" => args.w = value.parse()?,
            "s" => args.s = value.parse()?,
            "w1" => args.w1 = value.parse()?,
            "w2" => args.w2 = value.parse()?,
            "t" => args.t = value.parse()?,
            "stride" => args.stride = value.parse()?,
            "window" => args.window = value.parse()?,
            "subseq" => args.subseq = value.parse()?,
            "d" => args.d = value.parse()?,
            "repeats" => args.repeats = value.parse()?,
            "masks" => args.masks = value,
            other => bail!("未知参数：{other}"),
        }
        i += 2;
    }
    if args.input.is_empty() || args.output.is_empty() {
        bail!("必须提供 --input 与 --output");
    }
    if args.stride == 0 {
        bail!("--stride 必须 >= 1");
    }
    for mask in args.masks.split(',') {
        let mask = mask.trim();
        if mask.is_empty() || !mask.chars().all(|c| c == '0' || c == '1') {
            bail!("--masks 必须为逗号分隔的 0/1 串：{mask}");
        }
        if !mask.contains('1') {
            bail!("--masks 每个形状至少含一个 1：{mask}");
        }
    }
    Ok(args)
}

fn main() -> Result<()> {
    let args = parse_args()?;
    let records = read_fasta(&args.input)?;
    let started = Instant::now();
    let out = File::create(&args.output)?;
    let mut writer = BufWriter::new(out);
    let mut total_keys = 0usize;
    let mut total_residues = 0usize;

    for (id, seq) in &records {
        let cleaned = clean_protein(seq, SecPolicy::C);
        total_residues += cleaned.len();
        let mut rows: Vec<(u128, i64)> = Vec::new();
        match args.method.as_str() {
            "kmer" => {
                for (pos, key) in kmers(&cleaned, args.k)? {
                    if pos % args.stride == 0 {
                        rows.push((key, pos as i64));
                    }
                }
            }
            "minimizer" => {
                for key in minimizers(&cleaned, args.k, args.w)? {
                    rows.push((key, -1));
                }
            }
            "syncmer" => {
                for key in syncmers_with_mode(&cleaned, args.k, args.s, SyncmerMode::Closed)? {
                    rows.push((key, -1));
                }
            }
            "strobemer" => {
                for key in strobemers(&cleaned, args.k, args.w1, args.w2, args.t)? {
                    rows.push((key, -1));
                }
            }
            "subseqhash2" => {
                let params = SubseqHash2Params {
                    window: args.window,
                    subseq: args.subseq,
                    d: args.d,
                    repeats: args.repeats,
                    table_seed: 42,
                    chunk_size: 500,
                };
                for key in subseqhash2_keys(seq, &params)? {
                    rows.push((key, -1));
                }
            }
            "spaced" => {
                // 多形状间隔种子并集：按形状序号混盐防跨形状碰撞（i=0 与旧版单形状一致）
                for (i, mask) in args.masks.split(',').enumerate() {
                    let mask = mask.trim();
                    let weight = mask.chars().filter(|&c| c == '1').count();
                    let salt = (i as u128).wrapping_mul(0x9E3779B97F4A7C15);
                    for (pos, key) in spaced_kmers(&cleaned, mask, weight)? {
                        rows.push((key ^ salt, pos as i64));
                    }
                }
            }
            other => bail!("未知方法：{other}（本包仅含六类基线：kmer/spaced/minimizer/syncmer/strobemer/subseqhash2）"),
        }

        // 同一序列内按键去重（保留首次位置），保证集合语义。
        let mut seen = BTreeSet::new();
        let mut kept = 0usize;
        for (key, pos) in rows {
            if seen.insert(key) {
                writeln!(writer, "{id}\t{key:032x}\t{pos}")?;
                kept += 1;
            }
        }
        total_keys += kept;
    }
    writer.flush()?;
    let elapsed = started.elapsed().as_secs_f64() * 1000.0;
    eprintln!(
        "method={} seqs={} residues={} keys={} elapsed_ms={:.3}",
        args.method,
        records.len(),
        total_residues,
        total_keys,
        elapsed
    );
    Ok(())
}
