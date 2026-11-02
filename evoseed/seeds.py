"""蛋白质种子方法库：被阶段3分析与阶段5原型共用。

实现（均为蛋白质一级序列可用，不依赖结构）：
  - 氨基酸编码与物理化学替换等价类（进化容错）
  - 连续 k-mer、spaced k-mer（通配掩码）
  - minimizer（窗口最小哈希）与**标准 syncmer**（Edgar 2021：k-mer 内最小 s-mer 位于
    固定位置；open = 首端或末端，closed = 首端）
  - strobemer（anchor 配对，跨窗组合）
哈希为 128-bit（两个 splitmix64 组合），碰撞可忽略；跨进程稳定（不依赖 PYTHONHASHSEED）。
"""
from __future__ import annotations

from typing import Iterable

# 标准 20 氨基酸 + X 兜底
AA20 = "ACDEFGHIKLMNPQRSTVWY"
_AA_INDEX = {a: i for i, a in enumerate(AA20)}
# 等价类小写字母 a-f 的独立编码（接在 21..26），避免被 encode 误当 X
_CLASS_INDEX = {chr(ord("a") + i): 21 + i for i in range(6)}

# 残基策略（显式、可测）：
#  - 空格/制表/换行/'-'/'.' 视为空位，移除
#  - U（硒代半胱氨酸）默认归一为 C（Sec/Cys 氧化还原同源），可用 sec_policy="X" 改为未知
#  - '*'（终止）映射为 X
#  - 其余非 20 标准残基（X/B/Z/J/O 等）映射为 X
SEC_POLICIES = ("C", "X")


def clean_protein(seq: str, sec_policy: str = "C") -> str:
    """按显式规则清洗蛋白序列（见模块头部的残基策略）。"""
    if sec_policy not in SEC_POLICIES:
        raise ValueError(f"sec_policy must be one of {SEC_POLICIES}")
    out = []
    for ch in (seq or "").upper():
        if ch in (" ", "\t", "\r", "\n", "-", "."):
            continue
        if ch == "U":
            ch = "C" if sec_policy == "C" else "X"
        if ch == "*":
            ch = "X"
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


# 物理化学替换等价类（BLOSUM 启发，6 类；**两两互斥**，每个残基只属一类）
# 注：G（甘氨酸）归入 small_polar（b），不再同时属于 hydrophobic。
AA_CLASS = {
    "A": "a", "V": "a", "L": "a", "I": "a", "M": "a", "C": "a",
    "S": "b", "T": "b", "N": "b", "P": "b", "G": "b",
    "D": "c", "E": "c", "Q": "c",
    "K": "d", "R": "d", "H": "d",
    "F": "e", "Y": "e", "W": "e",
    "X": "f", "B": "f", "Z": "f", "J": "f", "O": "f", "U": "f",
}
SUB_CLASSES = {
    "hydrophobic": set("AVLIMC"),
    "small_polar": set("STNPG"),
    "acidic": set("DEQ"),
    "basic": set("KRH"),
    "aromatic": set("FYW"),
    "special": set("XBZJOU"),
}
_CLASS_OF = dict(AA_CLASS)

# 完整性自检：20 种标准氨基酸必须各自恰好属于一个类
assert set(AA_CLASS) >= set(AA20)


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
    """返回 128-bit 哈希（两个 splitmix64 拼接），使碰撞概率可忽略。

    先对残基编码做 64-bit 滚动词哈希，再用不同常量生成两个 64-bit 分量。
    """
    h = salt & 0xFFFFFFFFFFFFFFFF
    for c in codes:
        h = ((h << 5) + h + c + 1) & 0xFFFFFFFFFFFFFFFF
    lo = _splitmix64(h ^ 0x9E3779B97F4A7C15)
    hi = _splitmix64(h ^ 0xD1B54A32D192ED03)
    return (hi << 64) | lo


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


def syncmers(seq: str, k: int, s: int, mode: str = "open") -> set[int]:
    """标准 syncmer（Edgar 2021）。

    对每个长度 k 的 k-mer，计算其全部 (k-s+1) 个 s-mer 哈希；取最小哈希的位置 t
    （平局取最左）。若 t 落在规定位置，则保留该 k-mer：
      - mode="open"  : t == 0 或 t == k-s
      - mode="closed": t == 0
    返回被保留 k-mer 的 128-bit 哈希集合。
    """
    if k <= 0 or s <= 0 or s >= k:
        raise ValueError("require 0 < s < k")
    if mode not in ("open", "closed"):
        raise ValueError("mode must be 'open' or 'closed'")
    codes = encode(seq)
    n = len(codes)
    if n < k:
        return set()
    out = set()
    for i in range(n - k + 1):
        smers = [hash_kmer(codes[i + j:i + j + s]) for j in range(0, k - s + 1)]
        mn = min(smers)
        t = smers.index(mn)          # 平局取最左
        keep = (t == 0) if mode == "closed" else (t == 0 or t == k - s)
        if keep:
            out.add(hash_kmer(codes[i:i + k]))
    return out


def _anchors(seq: str, k: int, w1: int) -> list[tuple[int, int]]:
    """取 anchor 候选：窗口 w1 内的 minimizer（去重且保持位置顺序）。

    注：此处用 "anchor" 而非 "strobe"，以与 strobemer 的定义区分。
    """
    km = kmers(seq, k)
    if not km:
        return []
    seen = set()
    anchors = []
    for i in range(len(km) - w1 + 1):
        best = min(km[i:i + w1], key=lambda t: t[1])
        if best[1] not in seen:
            seen.add(best[1])
            anchors.append(best)
    return anchors


# 兼容旧名
_strobes = _anchors


def strobemers(seq: str, k: int, w1: int, w2: int, t: int = 2) -> set[int]:
    """**简化 strobemer 基线（fixed-offset anchor pair）**。

    取 minimizer anchors，把第 i 个 anchor 与其后第 t 个（距离在 [w1,w2]）配对。
    这不是标准的 randstrobe/minstrobe（后者在窗口内按哈希/随机选择第二 anchor），
    仅作为“固定偏移”对照基线；论文与 README 均按此标注。
    """
    anchors = _anchors(seq, k, w1)
    out = set()
    m = len(anchors)
    for i in range(m):
        p1, h1 = anchors[i]
        for j in range(i + 1, min(m, i + 6)):
            p2, h2 = anchors[j]
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
        return syncmers(s, p["k"], p["s"], p.get("mode", "open"))
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
