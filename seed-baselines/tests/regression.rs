//! 回归 golden：与论文导出（tools/seed_dump.exe 生成）逐字节比对。
//! golden 文件由论文在用二进制对 tests/fixtures/mini.faa 以 toolB 主口径参数导出。

use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_seed-baselines");
const FAA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/mini.faa");

fn dump(method: &[&str]) -> String {
    let out = std::env::temp_dir().join(format!("sb_golden_{}.tsv", method[1]));
    let st = Command::new(BIN)
        .arg("--input")
        .arg(FAA)
        .arg("--output")
        .arg(&out)
        .args(method)
        .status()
        .unwrap();
    assert!(st.success());
    std::fs::read_to_string(&out).unwrap()
}

#[test]
fn kmer_blastp_parity() {
    assert_eq!(
        dump(&["--method", "kmer", "--k", "3", "--stride", "1"]),
        include_str!("fixtures/mini_kmer.tsv")
    );
}

#[test]
fn kmer_diamond_parity() {
    assert_eq!(
        dump(&["--method", "spaced", "--masks", "111101110111,111011010010111"]),
        include_str!("fixtures/mini_kmer_diamond.tsv")
    );
}

#[test]
fn minimizer_parity() {
    assert_eq!(
        dump(&["--method", "minimizer", "--k", "7", "--w", "10"]),
        include_str!("fixtures/mini_minimizer.tsv")
    );
}

#[test]
fn syncmer_parity() {
    assert_eq!(
        dump(&["--method", "syncmer", "--k", "7", "--s", "5"]),
        include_str!("fixtures/mini_syncmer.tsv")
    );
}

#[test]
fn strobemer_parity() {
    assert_eq!(
        dump(&["--method", "strobemer", "--k", "7", "--w1", "9", "--w2", "32", "--t", "2"]),
        include_str!("fixtures/mini_strobemer.tsv")
    );
}

#[test]
fn subseqhash2_parity() {
    assert_eq!(
        dump(&["--method", "subseqhash2", "--window", "25", "--subseq", "11",
               "--d", "5", "--repeats", "11"]),
        include_str!("fixtures/mini_subseqhash2.tsv")
    );
}
