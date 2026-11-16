use anyhow::{bail, Context, Result};
use std::collections::BTreeMap;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// 严格读取 FASTA：拒绝空 header、重复 ID、空记录和无序列输入。
/// BufRead 避免在序列映射之外再常驻一份完整文件字符串。
pub fn read_fasta<P: AsRef<Path>>(path: P) -> Result<BTreeMap<String, String>> {
    let path = path.as_ref();
    let file =
        fs::File::open(path).with_context(|| format!("无法读取 FASTA：{}", path.display()))?;
    parse_fasta_reader(BufReader::new(file))
        .with_context(|| format!("FASTA 格式错误：{}", path.display()))
}

pub fn parse_fasta(text: &str) -> Result<BTreeMap<String, String>> {
    parse_fasta_reader(text.as_bytes())
}

pub fn parse_fasta_reader<R: BufRead>(mut reader: R) -> Result<BTreeMap<String, String>> {
    let mut records = BTreeMap::new();
    let mut current: Option<String> = None;
    let mut buffer = String::new();
    let mut raw = String::new();
    let mut line_no = 0_usize;

    loop {
        raw.clear();
        if reader.read_line(&mut raw)? == 0 {
            break;
        }
        line_no += 1;
        let line = raw.trim();
        if let Some(header) = line.strip_prefix('>') {
            flush_record(&mut records, &mut current, &mut buffer, line_no)?;
            let id = header.split_whitespace().next().unwrap_or("");
            if id.is_empty() {
                bail!("第 {line_no} 行：空 header");
            }
            current = Some(id.to_string());
        } else if !line.is_empty() {
            if current.is_none() {
                bail!("第 {line_no} 行：序列出现在首个 header 之前");
            }
            buffer.push_str(line);
        }
    }
    flush_record(&mut records, &mut current, &mut buffer, line_no + 1)?;
    if records.is_empty() {
        bail!("未找到任何 FASTA 记录");
    }
    Ok(records)
}

fn flush_record(
    records: &mut BTreeMap<String, String>,
    current: &mut Option<String>,
    buffer: &mut String,
    line_no: usize,
) -> Result<()> {
    let Some(id) = current.take() else {
        return Ok(());
    };
    if buffer.trim().is_empty() {
        bail!("记录 {id} 在约第 {line_no} 行之前没有序列");
    }
    if records.contains_key(&id) {
        bail!("重复序列 ID：{id}");
    }
    records.insert(id, std::mem::take(buffer));
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_duplicate_empty_and_orphan_records() {
        assert!(parse_fasta(">a\nAC\n>a\nDE\n")
            .unwrap_err()
            .to_string()
            .contains("重复"));
        assert!(parse_fasta(">\nAC\n")
            .unwrap_err()
            .to_string()
            .contains("空 header"));
        assert!(parse_fasta("AC\n").is_err());
        assert!(parse_fasta(">a\n\n").is_err());
        assert!(parse_fasta("").is_err());
    }
}
