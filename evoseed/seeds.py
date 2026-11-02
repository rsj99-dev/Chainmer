"""蛋白质种子方法库：被阶段3分析与阶段5原型共用。

实现（均为蛋白质一级序列可用，不依赖结构）：
  - 氨基酸编码与物理化学替换等价类（进化容错）
  - 连续 k-mer、spaced k-mer（通配掩码）
  - minimizer / syncmer（窗口内哈希选择）
  - strobemer（strobes 配对，跨窗组合）
哈希确定、跨进程稳定（不依赖 PYTHONHASHSEED）。
"""
from __future__ import annotations

from typing import Iterable

# 标准 20 氨基酸 + X 兜底
AA20 = "ACDEFGHIKLMNPQRSTVWY"
_AA_INDEX = {a: i for i, a in enumerate(AA20)}
_U = 20          # Sec(U) -> 视作 C 的氧化还原同源（编码上单独标记）
# 等价类小写字母 a-f 的独立编码（接在 21..26），避免被 encode 误当 X
_CLASS_INDEX = {chr(ord("a") + i): 21 + i for i in range(6)}


def clean_protein(seq: str) -> str:
    seq = seq.upper().strip()
    seq = seq.replace("U", "C")        # Sec 归一到 Cys
    out = []
    for ch in seq:
        out.append(ch if ch in _AA_INDEX else "X")
    return "".join(out)


def encode(seq: str) -> list[int]:
    codes = []
    for ch in seq:
        if ch in _AA_INDEX:
            codes.append(_AA_INDEX[ch])
        elif ch in _CLASS_INDEX:
            codes.append(_CLASS_INDEX[ch])
        else:
            codes.append(20)     # X 兼底
    return codes


# 物理化学替换等价类（BLOSUM 启发，6 类）——用于"进化容错型"种子
SUB_CLASSES = {
    "hydrophobic": set("AVLIMCG"),
    "small_polar": set("STNPG"),
    "acidic": set("DEQ"),
    "basic": set("KRH"),
    "aromatic": set("FYW"),
    "special": set("X"),
}
_CLASS_OF = {}
for _name, _s in SUB_CLASSES.items():
    for _a in _s:
        _CLASS_OF[_a] = {
            "hydrophobic": "a", "small_polar": "b", "acidic": "c",
            "basic": "d", "aromatic": "e", "special": "f",
        }[_name]


def to_class_string(seq: str) -> str:
    """把氨基酸序列映射为等价类串（替换容错：同类互替换不变）。"""
    return "".join(_CLASS_OF.get(ch, "f") for ch in seq)


# --------------------------------------------------------------------------- #
# 稳定哈希
# --------------------------------------------------------------------------- #
def _splitmix64(x: int) -> int:
    x = (x + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    return x ^ (x >> 31)


def hash_kmer(codes: Iterable[int], salt: int = 0xABCDEF) -> int:
    h = salt & 0xFFFFFFFFFFFFFFFF
    for c in codes:
        h = ((h << 5) + h + c + 1) & 0xFFFFFFFFFFFFFFFF
    return _splitmix64(h)


# --------------------------------------------------------------------------- #
# 种子提取
# --------------------------------------------------------------------------- #
def kmers(seq: str, k: int) -> list[tuple[int, int]]:
    """返回 [(pos, hash)]，连续 k-mer。"""
    codes = encode(seq)
    n = len(codes)
    if k <= 0 or n < k:
        return []
    out = []
    for i in range(n - k + 1):
        out.append((i, hash_kmer(codes[i:i + k])))
    return out


def spaced_kmers(seq: str, pattern: str, k: int) -> list[tuple[int, int]]:
    """spaced seed：pattern 为 '1'/'0'（或 '*' 通配）字符串，只对 '1' 位编码后哈希。

    k 指 '1' 的个数（有效位）。通配位（'0'）允许任意氨基酸，提高容错。
    """
    ones = [i for i, ch in enumerate(pattern) if ch == "1"]
    span = len(pattern)
    codes = encode(seq)
    n = len(codes)
    out = []
    for i in range(n - span + 1):
        sel = [codes[i + off] for off in ones]
        out.append((i, hash_kmer(sel)))
    return out


def minimizers(seq: str, k: int, w: int) -> set[int]:
    """窗口 w 上取哈希最小的 k-mer 作为 minimizer（去重集合）。"""
    km = kmers(seq, k)
    if not km:
        return set()
    out = set()
    for i in range(len(km) - w + 1):
        window = km[i:i + w]
        best = min(window, key=lambda t: t[1])
        out.add(best[1])
    return out


def syncmers(seq: str, k: int, w: int, s: int) -> set[int]:
    """syncmer：k-mer 是窗口内哈希最小且位于窗口内固定偏移 s 处，则选为标记。

    相比 minimizer，syncmer 由"内容"而非"位置"决定，保留更稳定。
    """
    km = kmers(seq, k)
    if not km or s >= w:
        return set()
    out = set()
    for i in range(len(km) - w + 1):
        window = km[i:i + w]
        # 最小哈希（平局取最左）
        minhash = min(h for _p, h in window)
        # 该窗口左起偏移 s 的 k-mer
        off_pos, off_hash = window[s]
        if off_hash == minhash and all(
                not (h < off_hash) for _p, h in window[:s]):
            out.add(off_hash)
    return out


def _strobes(seq: str, k: int, w1: int) -> list[tuple[int, int]]:
    """取 strobe 候选：窗口 w1 内的 minimizer（去重且保持位置顺序）。"""
    km = kmers(seq, k)
    if not km:
        return []
    seen = set()
    strobes = []
    for i in range(len(km) - w1 + 1):
        best = min(km[i:i + w1], key=lambda t: t[1])
        if best[1] not in seen:
            seen.add(best[1])
            strobes.append(best)
    return strobes


def strobemers(seq: str, k: int, w1: int, w2: int, t: int = 2) -> set[int]:
    """strobemer：把第 i 个 strobe 与其后第 t 个、距离落在 [w1, w2] 的 strobe 组合哈希。"""
    strobes = _strobes(seq, k, w1)
    out = set()
    m = len(strobes)
    for i in range(m):
        p1, h1 = strobes[i]
        for j in range(i + 1, min(m, i + 6)):
            p2, h2 = strobes[j]
            d = p2 - p1
            if d < w1:
                continue
            if d > w2:
                break
            if (j - i) == t or t <= 0:
                out.add(hash_kmer([h1, h2], salt=0x51ED))
    return out


# --------------------------------------------------------------------------- #
# 统一接口：按规格抽取种子集合
# --------------------------------------------------------------------------- #
SEED_METHODS = ("kmer", "spaced", "minimizer", "syncmer", "strobemer",
                "ckmer", "cspaced")   # c* = 替换等价类版本（进化容错）


def extract_seeds(seq: str, method: str, **p) -> set[int]:
    s = clean_protein(seq)
    if method == "ckmer":
        s = to_class_string(s)
    elif method == "cspaced":
        s = to_class_string(s)
    if method in ("kmer", "ckmer"):
        return {h for _p, h in kmers(s, p["k"])}
    if method in ("spaced", "cspaced"):
        return {h for _p, h in spaced_kmers(s, p["pattern"], p.get("k", 0))}
    if method == "minimizer":
        return minimizers(s, p["k"], p["w"])
    if method == "syncmer":
        return syncmers(s, p["k"], p["w"], p["s"])
    if method == "strobemer":
        return strobemers(s, p["k"], p["w1"], p["w2"], p.get("t", 2))
    raise ValueError(f"unknown method {method}")


def seed_positions(seq: str, method: str, **p) -> list[tuple[int, int]]:
    """返回 (pos, hash) 列表，用于需要位置信息的组合种子。"""
    s = clean_protein(seq)
    if method in ("ckmer", "cspaced"):
        s = to_class_string(s)
    if method in ("kmer", "ckmer"):
        return kmers(s, p["k"])
    if method in ("spaced", "cspaced"):
        return spaced_kmers(s, p["pattern"], p.get("k", 0))
    return []
