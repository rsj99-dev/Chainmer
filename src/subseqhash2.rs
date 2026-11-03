//! SubseqHash2 的蛋白字母表标量适配。
//!
//! 按 Shao-Group 官方 C++ 实现的 ABCk orders、前向/反向 DP 与 pivot 合并移植。
//! 官方实现面向四字符核酸并含 SIMD；本模块面向 20 种标准氨基酸 + X，不宣称
//! 与官方 SIMD 二进制具有相同运行时间。X 类似官方 N，包含 X 的窗口不产种子。

use crate::seeds::{clean_protein, encode, hash_kmer, SecPolicy, AMBIGUOUS};
use anyhow::{bail, Result};
use std::collections::BTreeSet;

const INF: i64 = 1_i64 << 62;
const PROTEIN_SIGMA: usize = 21;
const SUBSEQ_SALT: u64 = 0x5EED2;

#[derive(Clone, Debug)]
pub struct SubseqHash2Params {
    pub window: usize,
    pub subseq: usize,
    pub d: usize,
    pub repeats: usize,
    pub table_seed: u64,
    pub chunk_size: usize,
}

impl Default for SubseqHash2Params {
    fn default() -> Self {
        // 固定、未在 OrthoChain 测试标签上调优的蛋白适配参数。
        Self {
            window: 25,
            subseq: 11,
            d: 5,
            repeats: 11,
            table_seed: 42,
            chunk_size: 500,
        }
    }
}

impl SubseqHash2Params {
    pub fn validate(&self) -> Result<()> {
        if self.window == 0 || self.subseq < 2 || self.window < self.subseq {
            bail!("SubseqHash2 要求 window >= subseq >= 2");
        }
        if self.d == 0 || self.d > 31 {
            bail!("SubseqHash2 要求 1 <= d <= 31");
        }
        if self.repeats == 0 || self.repeats > self.subseq {
            bail!("SubseqHash2 要求 1 <= repeats <= subseq");
        }
        if self.chunk_size < self.window {
            bail!("chunk_size 必须 >= window");
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Default)]
struct DpCell {
    f_min: i64,
    f_max: i64,
    g_min: usize,
    g_max: usize,
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn range(&mut self, low: i64, high: i64) -> i64 {
        low + (self.next() % (high - low + 1) as u64) as i64
    }
    fn sign(&mut self) -> i8 {
        if self.next() & 1 == 0 {
            -1
        } else {
            1
        }
    }
}

struct AbcTables {
    d: usize,
    sigma: usize,
    a: Vec<i64>,
    rev_a: Vec<i64>,
    b1: Vec<i8>,
    b2: Vec<i8>,
    rev_b1: Vec<i8>,
    rev_b2: Vec<i8>,
    c1: Vec<usize>,
    c2: Vec<usize>,
    c3: Vec<usize>,
    a3: Vec<i64>,
    combine1: Vec<i8>,
    combine2: Vec<i8>,
    combine3: Vec<i8>,
    valid: Vec<bool>,
}

impl AbcTables {
    fn new(params: &SubseqHash2Params) -> Self {
        let (k, d, sigma) = (params.subseq, params.d, PROTEIN_SIGMA);
        let mut rng = Rng(params.table_seed);
        let size3 = k * sigma * d;
        let size2 = k * sigma;
        let mut a = vec![0; size3];
        let mut rev_a = vec![0; size3];
        let mut b1 = vec![0; size3];
        let mut b2 = vec![0; size3];
        let mut rev_b1 = vec![0; size3];
        let mut rev_b2 = vec![0; size3];
        let mut c1 = vec![0; size2];
        let mut c2 = vec![0; size2];
        let mut c3 = vec![0; size2];
        let mut a3 = vec![0; size2];
        let mut combine1 = vec![0; size2];
        let mut combine2 = vec![0; size2];
        let mut combine3 = vec![0; size2];
        for i in 0..k {
            for j in 0..sigma {
                for q in 0..d {
                    let idx = (i * sigma + j) * d + q;
                    a[idx] = rng.range(10, 100);
                    b1[idx] = rng.sign();
                    b2[idx] = rng.sign();
                }
            }
        }
        for i in 0..k {
            for j in 0..sigma {
                c1[i * sigma + j] = (rng.next() % d as u64) as usize;
            }
        }
        for idx in 0..size3 {
            rev_a[idx] = rng.range(10, 100);
            rev_b1[idx] = rng.sign();
            rev_b2[idx] = rng.sign();
        }
        for i in 0..k {
            for j in 0..sigma {
                c2[i * sigma + j] = (rng.next() % d as u64) as usize;
            }
        }
        for idx in 0..size2 {
            combine1[idx] = rng.sign();
            combine2[idx] = rng.sign();
            a3[idx] = rng.range(10, 100);
            combine3[idx] = rng.sign();
            c3[idx] = (rng.next() % d as u64) as usize;
        }
        let mut order: Vec<usize> = (0..k).collect();
        for i in (1..k).rev() {
            order.swap(i, rng.next() as usize % (i + 1));
        }
        let mut valid = vec![false; k];
        for pivot in order.into_iter().take(params.repeats) {
            valid[pivot] = true;
        }
        Self {
            d,
            sigma,
            a,
            rev_a,
            b1,
            b2,
            rev_b1,
            rev_b2,
            c1,
            c2,
            c3,
            a3,
            combine1,
            combine2,
            combine3,
            valid,
        }
    }

    fn i3(&self, i: usize, nt: usize, q: usize) -> usize {
        (i * self.sigma + nt) * self.d + q
    }
    fn i2(&self, i: usize, nt: usize) -> usize {
        i * self.sigma + nt
    }
}

struct Layout {
    d: usize,
    dim1: usize,
    dim2: usize,
}
impl Layout {
    fn new(n: usize, k: usize, d: usize) -> Self {
        let dim2 = (k + 1) * d;
        let dim1 = (n + 1) * dim2;
        Self { d, dim1, dim2 }
    }
    fn dp(&self, st: usize, i: usize, j: usize, q: usize) -> usize {
        st * self.dim1 + i * self.dim2 + j * self.d + q
    }
    fn h(&self, i: usize, j: usize, q: usize) -> usize {
        i * self.dim2 + j * self.d + q
    }
}

fn forward(
    codes: &[u8],
    start: usize,
    end: usize,
    p: &SubseqHash2Params,
    t: &AbcTables,
) -> Vec<DpCell> {
    let len = end - start;
    let lay = Layout::new(p.window, p.subseq, p.d);
    let del = p.window - p.subseq;
    let mut dp = vec![DpCell::default(); lay.dim1 * len];
    let mut h = vec![0usize; lay.dim1];
    for st in 0..len {
        let en = (st + p.window - 1).min(len - 1);
        let interval = en - st + 1;
        for i in 0..=interval {
            let idx = lay.dp(st, i, 0, 0);
            dp[idx].f_min = 0;
            dp[idx].f_max = 0;
            h[lay.h(i, 0, 0)] = st + 1;
            for q in 1..p.d {
                let x = lay.dp(st, i, 0, q);
                dp[x].f_min = INF;
                dp[x].f_max = -INF;
            }
        }
        let nt0 = codes[start + st] as usize;
        let d1 = t.c1[t.i2(0, nt0)];
        for q in 0..p.d {
            let x = lay.dp(st, 1, 1, q);
            dp[x].f_min = INF;
            dp[x].f_max = -INF;
        }
        if nt0 as u8 == AMBIGUOUS {
            continue;
        }
        let idx = lay.dp(st, 1, 1, d1);
        let z = t.i3(0, nt0, d1);
        dp[idx].f_min = t.b2[z] as i64 * t.a[z];
        dp[idx].f_max = dp[idx].f_min;
        h[lay.h(1, 1, d1)] = st + 1;
        for i in 2..=interval {
            if codes[start + st + i - 1] == AMBIGUOUS {
                break;
            }
            let minj = 1.max(i.saturating_sub(del));
            let maxj = i.min(p.subseq - 1);
            for j in minj..=maxj {
                let nt = codes[start + st + i - 1] as usize;
                let v = t.c1[t.i2(j - 1, nt)];
                for q in 0..p.d {
                    let z = (q + v) % p.d;
                    let idx = lay.dp(st, i, j, z);
                    if h[lay.h(i, j, z)] != st + 1 {
                        dp[idx].f_min = INF;
                        dp[idx].f_max = -INF;
                    }
                    let prev = lay.dp(st, i - 1, j, z);
                    if h[lay.h(i - 1, j, z)] == st + 1 {
                        if dp[prev].f_min < dp[idx].f_min {
                            dp[idx].f_min = dp[prev].f_min;
                            dp[idx].g_min = dp[prev].g_min;
                        }
                        if dp[prev].f_max > dp[idx].f_max {
                            dp[idx].f_max = dp[prev].f_max;
                            dp[idx].g_max = dp[prev].g_max;
                        }
                        h[lay.h(i, j, z)] = st + 1;
                    }
                    let from = lay.dp(st, i - 1, j - 1, q);
                    if h[lay.h(i - 1, j - 1, q)] == st + 1 {
                        let tab = t.i3(j - 1, nt, z);
                        let av = t.a[tab] * (t.b2[tab] as i64);
                        if t.b1[tab] == -1 {
                            let v1 = -dp[from].f_min + av;
                            let v2 = -dp[from].f_max + av;
                            if v1 > dp[idx].f_max {
                                dp[idx].f_max = v1;
                                dp[idx].g_max = i - 1;
                            }
                            if v2 < dp[idx].f_min {
                                dp[idx].f_min = v2;
                                dp[idx].g_min = i - 1;
                            }
                        } else {
                            let v1 = dp[from].f_min + av;
                            let v2 = dp[from].f_max + av;
                            if v1 < dp[idx].f_min {
                                dp[idx].f_min = v1;
                                dp[idx].g_min = i - 1;
                            }
                            if v2 > dp[idx].f_max {
                                dp[idx].f_max = v2;
                                dp[idx].g_max = i - 1;
                            }
                        }
                        h[lay.h(i, j, z)] = st + 1;
                    }
                }
            }
        }
    }
    dp
}

fn reverse(
    codes: &[u8],
    start: usize,
    end: usize,
    p: &SubseqHash2Params,
    t: &AbcTables,
) -> Vec<DpCell> {
    let len = end - start;
    let lay = Layout::new(p.window, p.subseq, p.d);
    let del = p.window - p.subseq;
    let mut dp = vec![DpCell::default(); lay.dim1 * len];
    let mut h = vec![0usize; lay.dim1];
    for st in (0..len).rev() {
        let en = st.saturating_sub(p.window - 1);
        let interval = st - en + 1;
        for i in 0..=interval {
            let idx = lay.dp(st, i, 0, 0);
            dp[idx].f_min = 0;
            dp[idx].f_max = 0;
            h[lay.h(i, 0, 0)] = st + 1;
            for q in 1..p.d {
                let x = lay.dp(st, i, 0, q);
                dp[x].f_min = INF;
                dp[x].f_max = -INF;
            }
        }
        let nt0 = codes[start + st] as usize;
        let d1 = t.c2[t.i2(0, nt0)];
        for q in 0..p.d {
            let x = lay.dp(st, 1, 1, q);
            dp[x].f_min = INF;
            dp[x].f_max = -INF;
        }
        if nt0 as u8 == AMBIGUOUS {
            continue;
        }
        let idx = lay.dp(st, 1, 1, d1);
        let z = t.i3(0, nt0, d1);
        dp[idx].f_min = t.rev_b2[z] as i64 * t.rev_a[z];
        dp[idx].f_max = dp[idx].f_min;
        h[lay.h(1, 1, d1)] = st + 1;
        for i in 2..=interval {
            if codes[start + st + 1 - i] == AMBIGUOUS {
                break;
            }
            let minj = 1.max(i.saturating_sub(del));
            let maxj = i.min(p.subseq - 1);
            for j in minj..=maxj {
                let nt = codes[start + st + 1 - i] as usize;
                let v = t.c2[t.i2(j - 1, nt)];
                for q in 0..p.d {
                    let z = (q + v) % p.d;
                    let idx = lay.dp(st, i, j, z);
                    if h[lay.h(i, j, z)] != st + 1 {
                        dp[idx].f_min = INF;
                        dp[idx].f_max = -INF;
                    }
                    let prev = lay.dp(st, i - 1, j, z);
                    if h[lay.h(i - 1, j, z)] == st + 1 {
                        if dp[prev].f_min < dp[idx].f_min {
                            dp[idx].f_min = dp[prev].f_min;
                            dp[idx].g_min = dp[prev].g_min;
                        }
                        if dp[prev].f_max > dp[idx].f_max {
                            dp[idx].f_max = dp[prev].f_max;
                            dp[idx].g_max = dp[prev].g_max;
                        }
                        h[lay.h(i, j, z)] = st + 1;
                    }
                    let from = lay.dp(st, i - 1, j - 1, q);
                    if h[lay.h(i - 1, j - 1, q)] == st + 1 {
                        let tab = t.i3(j - 1, nt, z);
                        let av = t.rev_a[tab] * (t.rev_b2[tab] as i64);
                        if t.rev_b1[tab] == -1 {
                            let v1 = -dp[from].f_min + av;
                            let v2 = -dp[from].f_max + av;
                            if v1 > dp[idx].f_max {
                                dp[idx].f_max = v1;
                                dp[idx].g_max = i - 1;
                            }
                            if v2 < dp[idx].f_min {
                                dp[idx].f_min = v2;
                                dp[idx].g_min = i - 1;
                            }
                        } else {
                            let v1 = dp[from].f_min + av;
                            let v2 = dp[from].f_max + av;
                            if v1 < dp[idx].f_min {
                                dp[idx].f_min = v1;
                                dp[idx].g_min = i - 1;
                            }
                            if v2 > dp[idx].f_max {
                                dp[idx].f_max = v2;
                                dp[idx].g_max = i - 1;
                            }
                        }
                        h[lay.h(i, j, z)] = st + 1;
                    }
                }
            }
        }
    }
    dp
}

fn better(v: i64, target: usize, best_v: i64, best_target: usize) -> bool {
    (v > best_v && target == best_target) || target < best_target
}

struct ChunkDp<'a> {
    codes: &'a [u8],
    start: usize,
    end: usize,
    dp: &'a [DpCell],
    rev: &'a [DpCell],
}

fn first_finite(
    dp: &[DpCell],
    lay: &Layout,
    st: usize,
    i: usize,
    j: usize,
    initial_q: usize,
    d: usize,
) -> Option<(usize, usize)> {
    for offset in 0..d {
        let q = (initial_q + offset) % d;
        let idx = lay.dp(st, i, j, q);
        if dp[idx].f_min < INF {
            return Some((q, idx));
        }
    }
    None
}

fn combine(
    chunk: ChunkDp<'_>,
    p: &SubseqHash2Params,
    t: &AbcTables,
    out: &mut BTreeSet<u128>,
) -> Result<()> {
    let ChunkDp {
        codes,
        start,
        end,
        dp,
        rev,
    } = chunk;
    let len = end - start;
    let lay = Layout::new(p.window, p.subseq, p.d);
    let del = p.window - p.subseq;
    for st in 0..=len - p.window {
        if codes[start + st..start + st + p.window].contains(&AMBIGUOUS) {
            continue;
        }
        let mut best_score = vec![-INF; p.subseq];
        let mut best_pos = vec![None; p.subseq];
        let mut best_q = vec![0usize; p.subseq];
        let mut best_psi = vec![p.d; p.subseq];
        for i in 0..p.window {
            let nt = codes[start + st + i] as usize;
            if p.window - i >= p.subseq && t.valid[0] {
                let d3 = t.c3[t.i2(0, nt)];
                let Some((d1, idx)) = first_finite(
                    dp,
                    &lay,
                    st + i + 1,
                    p.window - i - 1,
                    p.subseq - 1,
                    (p.d - d3) % p.d,
                    p.d,
                ) else {
                    continue;
                };
                let mut v = t.combine3[t.i2(0, nt)] as i64 * t.a3[t.i2(0, nt)];
                if t.combine1[t.i2(0, nt)] == 1 {
                    v += dp[idx].f_max;
                } else {
                    v -= dp[idx].f_min;
                }
                let psi = (d1 + d3) % p.d;
                if better(v, psi, best_score[0], best_psi[0]) {
                    best_score[0] = v;
                    best_pos[0] = Some(i);
                    best_q[0] = d1;
                    best_psi[0] = psi;
                }
            }
            if i >= p.subseq - 1 && t.valid[p.subseq - 1] {
                let d3 = t.c3[t.i2(p.subseq - 1, nt)];
                let Some((d1, idx)) = first_finite(
                    rev,
                    &lay,
                    st + i - 1,
                    i,
                    p.subseq - 1,
                    (p.d - d3) % p.d,
                    p.d,
                ) else {
                    continue;
                };
                let mut v =
                    t.combine3[t.i2(p.subseq - 1, nt)] as i64 * t.a3[t.i2(p.subseq - 1, nt)];
                if t.combine2[t.i2(p.subseq - 1, nt)] == 1 {
                    v += rev[idx].f_max;
                } else {
                    v -= rev[idx].f_min;
                }
                let psi = (d1 + d3) % p.d;
                let j = p.subseq - 1;
                if better(v, psi, best_score[j], best_psi[j]) {
                    best_score[j] = v;
                    best_pos[j] = Some(i);
                    best_q[j] = 0;
                    best_psi[j] = psi;
                }
            }
            let minj = 1.max(i.saturating_sub(del));
            let maxj = i.min(p.subseq.saturating_sub(2));
            if minj > maxj {
                continue;
            }
            for j in minj..=maxj {
                if !t.valid[j] {
                    continue;
                }
                let fbase = lay.dp(st + i + 1, p.window - i - 1, p.subseq - j - 1, 0);
                let rbase = lay.dp(st + i - 1, i, j, 0);
                let mut left_mask = 0u64;
                for q in 0..p.d {
                    if rev[rbase + q].f_min < INF {
                        left_mask |= 1u64 << q;
                    }
                }
                let mut total_mask = 0u64;
                for q in 0..p.d {
                    if dp[fbase + q].f_min < INF {
                        total_mask |= left_mask << q;
                    }
                }
                let d3 = t.c3[t.i2(j, nt)];
                let mut d1 = (p.d - d3) % p.d;
                let mut target = 0usize;
                while target < p.d {
                    if ((total_mask >> d1) & 1) != 0 || ((total_mask >> (d1 + p.d)) & 1) != 0 {
                        break;
                    }
                    d1 = (d1 + 1) % p.d;
                    target += 1;
                }
                if target == p.d {
                    continue;
                }
                for q in 0..p.d {
                    let fi = lay.dp(st + i + 1, p.window - i - 1, p.subseq - j - 1, q);
                    let ri = lay.dp(st + i - 1, i, j, (target + 2 * p.d - d3 - q) % p.d);
                    if dp[fi].f_min >= INF || rev[ri].f_min >= INF {
                        continue;
                    }
                    let mut v = t.combine3[t.i2(j, nt)] as i64 * t.a3[t.i2(j, nt)];
                    if t.combine1[t.i2(j, nt)] == 1 {
                        v += dp[fi].f_max;
                    } else {
                        v -= dp[fi].f_min;
                    }
                    if t.combine2[t.i2(j, nt)] == 1 {
                        v += rev[ri].f_max;
                    } else {
                        v -= rev[ri].f_min;
                    }
                    if better(v, target, best_score[j], best_psi[j]) {
                        best_score[j] = v;
                        best_pos[j] = Some(i);
                        best_q[j] = q;
                        best_psi[j] = target;
                    }
                }
            }
        }
        for j in 0..p.subseq {
            let Some(pivot) = best_pos[j] else {
                continue;
            };
            let nt = codes[start + st + pivot] as usize;
            let mut positions = Vec::with_capacity(p.subseq);
            let mut d2 = (best_psi[j] + 2 * p.d - best_q[j] - t.c3[t.i2(j, nt)]) % p.d;
            let mut x = pivot;
            let mut y = j;
            let mut z = t.combine2[t.i2(j, nt)] >= 0;
            while y > 0 {
                let idx = lay.dp(st + pivot - 1, x, y, d2);
                let g = if z { rev[idx].g_max } else { rev[idx].g_min };
                let pos = start + st + pivot - g - 1;
                positions.push(pos);
                let next = codes[pos] as usize;
                if t.rev_b1[t.i3(y - 1, next, d2)] == -1 {
                    z = !z;
                }
                d2 = (d2 + p.d - t.c2[t.i2(y - 1, next)]) % p.d;
                x = g;
                y -= 1;
            }
            positions.push(start + st + pivot);
            if j < p.subseq - 1 {
                let mut d1 = best_q[j];
                x = p.window - pivot - 1;
                y = p.subseq - j - 1;
                z = t.combine1[t.i2(j, nt)] >= 0;
                while y > 0 {
                    let idx = lay.dp(st + pivot + 1, x, y, d1);
                    let g = if z { dp[idx].g_max } else { dp[idx].g_min };
                    let pos = start + st + pivot + g + 1;
                    positions.push(pos);
                    let next = codes[pos] as usize;
                    if t.b1[t.i3(y - 1, next, d1)] == -1 {
                        z = !z;
                    }
                    d1 = (d1 + p.d - t.c1[t.i2(y - 1, next)]) % p.d;
                    x = g;
                    y -= 1;
                }
            }
            positions.sort_unstable();
            positions.dedup();
            if positions.len() != p.subseq {
                bail!("SubseqHash2 内部回溯错误：子序列长度不符");
            }
            let subseq: Vec<u64> = positions.into_iter().map(|pos| codes[pos] as u64).collect();
            out.insert(hash_kmer(&subseq, SUBSEQ_SALT));
        }
    }
    Ok(())
}

pub fn subseqhash2_keys(seq: &str, params: &SubseqHash2Params) -> Result<BTreeSet<u128>> {
    params.validate()?;
    let cleaned = clean_protein(seq, SecPolicy::C);
    let codes = encode(&cleaned);
    if codes.len() < params.window {
        return Ok(BTreeSet::new());
    }
    let tables = AbcTables::new(params);
    let mut out = BTreeSet::new();
    let mut start = 0usize;
    loop {
        let end = start.saturating_add(params.chunk_size).min(codes.len());
        let dp = forward(&codes, start, end, params, &tables);
        let rev = reverse(&codes, start, end, params, &tables);
        combine(
            ChunkDp {
                codes: &codes,
                start,
                end,
                dp: &dp,
                rev: &rev,
            },
            params,
            &tables,
            &mut out,
        )?;
        if end == codes.len() {
            break;
        }
        let next = end - params.window + 1;
        if next <= start {
            bail!("SubseqHash2 分块推进失败");
        }
        start = next;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn small() -> SubseqHash2Params {
        SubseqHash2Params {
            window: 8,
            subseq: 4,
            d: 3,
            repeats: 4,
            table_seed: 7,
            chunk_size: 16,
        }
    }

    fn naive_forward_score(
        codes: &[u8],
        start: usize,
        st: usize,
        state: (usize, usize, usize),
        p: &SubseqHash2Params,
        t: &AbcTables,
    ) -> Option<(i64, i64)> {
        let (i, j, q) = state;
        if j == 0 {
            return (q == 0).then_some((0, 0));
        }
        if i == 0 || j > i || codes[start + st..start + st + i].contains(&AMBIGUOUS) {
            return None;
        }
        let mut scores = Vec::new();
        for mask in 0usize..(1usize << i) {
            if mask.count_ones() as usize != j {
                continue;
            }
            let mut state = 0usize;
            let mut score = 0i64;
            let mut selected = 0usize;
            for pos in 0..i {
                if mask & (1usize << pos) == 0 {
                    continue;
                }
                let nt = codes[start + st + pos] as usize;
                let z = (state + t.c1[t.i2(selected, nt)]) % p.d;
                let tab = t.i3(selected, nt, z);
                score = t.b1[tab] as i64 * score + t.b2[tab] as i64 * t.a[tab];
                state = z;
                selected += 1;
            }
            if state == q {
                scores.push(score);
            }
        }
        let min = scores.iter().copied().min()?;
        let max = scores.iter().copied().max()?;
        Some((min, max))
    }

    fn naive_reverse_score(
        codes: &[u8],
        start: usize,
        st: usize,
        state: (usize, usize, usize),
        p: &SubseqHash2Params,
        t: &AbcTables,
    ) -> Option<(i64, i64)> {
        let (i, j, q) = state;
        if j == 0 {
            return (q == 0).then_some((0, 0));
        }
        if i == 0 || j > i || (0..i).any(|pos| codes[start + st - pos] == AMBIGUOUS) {
            return None;
        }
        let mut scores = Vec::new();
        for mask in 0usize..(1usize << i) {
            if mask.count_ones() as usize != j {
                continue;
            }
            let mut state = 0usize;
            let mut score = 0i64;
            let mut selected = 0usize;
            for pos in 0..i {
                if mask & (1usize << pos) == 0 {
                    continue;
                }
                let nt = codes[start + st - pos] as usize;
                let z = (state + t.c2[t.i2(selected, nt)]) % p.d;
                let tab = t.i3(selected, nt, z);
                score = t.rev_b1[tab] as i64 * score + t.rev_b2[tab] as i64 * t.rev_a[tab];
                state = z;
                selected += 1;
            }
            if state == q {
                scores.push(score);
            }
        }
        let min = scores.iter().copied().min()?;
        let max = scores.iter().copied().max()?;
        Some((min, max))
    }

    #[test]
    fn reverse_dp_matches_naive_subsequence_enumeration() {
        let p = SubseqHash2Params {
            window: 6,
            subseq: 4,
            d: 3,
            repeats: 4,
            table_seed: 19,
            chunk_size: 16,
        };
        let cleaned = clean_protein("ACDEFGHI", SecPolicy::C);
        let codes = encode(&cleaned);
        let tables = AbcTables::new(&p);
        let dp = reverse(&codes, 0, codes.len(), &p, &tables);
        let lay = Layout::new(p.window, p.subseq, p.d);
        for st in 0..codes.len() {
            let interval = p.window.min(st + 1);
            for i in 0..=interval {
                for j in 0..=i.min(p.subseq) {
                    let minj = 1.max(i.saturating_sub(p.window - p.subseq));
                    let maxj = i.min(p.subseq - 1);
                    let intended =
                        j == 0 || (i == 1 && j == 1) || (i >= 2 && j >= minj && j <= maxj);
                    if !intended {
                        continue;
                    }
                    for q in 0..p.d {
                        let cell = dp[lay.dp(st, i, j, q)];
                        match naive_reverse_score(&codes, 0, st, (i, j, q), &p, &tables) {
                            Some((min, max)) => {
                                assert_eq!(
                                    (cell.f_min, cell.f_max),
                                    (min, max),
                                    "反向 DP 数值不符：st={st} i={i} j={j} q={q}"
                                );
                            }
                            None => {
                                assert!(
                                    cell.f_min >= INF && cell.f_max <= -INF,
                                    "反向 DP 应不可达：st={st} i={i} j={j} q={q}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn forward_dp_matches_naive_subsequence_enumeration() {
        let p = SubseqHash2Params {
            window: 6,
            subseq: 4,
            d: 3,
            repeats: 4,
            table_seed: 19,
            chunk_size: 16,
        };
        let cleaned = clean_protein("ACDEFGHI", SecPolicy::C);
        let codes = encode(&cleaned);
        let tables = AbcTables::new(&p);
        let dp = forward(&codes, 0, codes.len(), &p, &tables);
        let lay = Layout::new(p.window, p.subseq, p.d);
        for st in 0..codes.len() {
            let interval = p.window.min(codes.len() - st);
            for i in 0..=interval {
                for j in 0..=i.min(p.subseq) {
                    let minj = 1.max(i.saturating_sub(p.window - p.subseq));
                    let maxj = i.min(p.subseq - 1);
                    let intended =
                        j == 0 || (i == 1 && j == 1) || (i >= 2 && j >= minj && j <= maxj);
                    if !intended {
                        continue;
                    }
                    for q in 0..p.d {
                        let cell = dp[lay.dp(st, i, j, q)];
                        match naive_forward_score(&codes, 0, st, (i, j, q), &p, &tables) {
                            Some((min, max)) => {
                                assert_eq!(
                                    (cell.f_min, cell.f_max),
                                    (min, max),
                                    "DP 数值不符：st={st} i={i} j={j} q={q}"
                                );
                            }
                            None => {
                                assert!(
                                    cell.f_min >= INF && cell.f_max <= -INF,
                                    "DP 应不可达：st={st} i={i} j={j} q={q}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn deterministic_bounded_and_short_sequences() {
        let p = small();
        let seq = "ACDEFGHIKLMNPQRSTVWY";
        let a = subseqhash2_keys(seq, &p).unwrap();
        let b = subseqhash2_keys(seq, &p).unwrap();
        assert_eq!(a, b);
        assert!(!a.is_empty());
        assert!(a.len() <= (seq.len() - p.window + 1) * p.repeats);
        assert!(subseqhash2_keys("ACD", &p).unwrap().is_empty());
    }

    #[test]
    fn ambiguous_windows_and_invalid_parameters() {
        let p = small();
        assert!(subseqhash2_keys("ACDEXXXXACDE", &p).unwrap().is_empty());
        let mut bad = small();
        bad.d = 0;
        assert!(subseqhash2_keys("ACDEFGHI", &bad).is_err());
        bad = small();
        bad.repeats = 99;
        assert!(subseqhash2_keys("ACDEFGHI", &bad).is_err());
        bad = small();
        bad.subseq = 1;
        assert!(subseqhash2_keys("ACDEFGHI", &bad).is_err());
    }

    #[test]
    fn chunking_does_not_change_keys() {
        let seq = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY";
        for chunk in [8, 9, 12, 17, 64] {
            let mut a = small();
            a.chunk_size = chunk;
            let mut b = small();
            b.chunk_size = 64;
            assert_eq!(
                subseqhash2_keys(seq, &a).unwrap(),
                subseqhash2_keys(seq, &b).unwrap()
            );
        }
    }

    #[test]
    fn minimal_window_and_single_modulo_state_are_supported() {
        let p = SubseqHash2Params {
            window: 2,
            subseq: 2,
            d: 1,
            repeats: 1,
            table_seed: 3,
            chunk_size: 2,
        };
        assert_eq!(subseqhash2_keys("AC", &p).unwrap().len(), 1);
        assert_eq!(subseqhash2_keys("ACDE", &p).unwrap().len(), 3);
    }
}
