"""EvoSeed 原型实现 —— 方案一 COSSEED。

种子 = 替换等价类空间中的"有序锚点对（anchor pair）"，编码：
    K = H( h(c_kmer_a), h(c_kmer_b), bucket(gap) )
其中 c_kmer 为类串上的 kc 长度子串；用类 minimizer 作为 anchor 以控制数量；
两 anchor 的间距落在 [gmin,gmax] 且按 delta 分桶，保证"离散保守位点组合 + 相对间距一致性"。
哈希为 128-bit（见 seedlib.hash_kmer），碰撞可忽略。

可消融开关：
  use_class : False 时退化为原始氨基酸（去掉进化容错）
  use_gap   : False 时不记录间距桶（退化为纯 k-mer 对，去掉位置一致性）
  w1        : minimizer 窗口，控制 anchor 密度
匹配分数 = 两条序列共享的锚点对键数（>=m 视为候选命中，m 为多重一致性阈值）。
"""
from __future__ import annotations

from . import seeds as seedlib


def _class_or_raw(seq: str, use_class: bool) -> str:
    s = seedlib.clean_protein(seq)
    return seedlib.to_class_string(s) if use_class else s


def anchors(seq: str, kc: int, w1: int, use_class: bool) -> list[tuple[int, int]]:
    """类(或原始)串上的 minimizer 作为 anchor：窗口 w1 内哈希最小的 c_kmer。

    返回按位置递增的 anchor 列表；并列时 min() 取窗口内首个（最左）。
    """
    sc = _class_or_raw(seq, use_class)
    codes = seedlib.encode(sc)
    n = len(codes)
    if n < kc:
        return []
    km = []
    for i in range(n - kc + 1):
        km.append((i, seedlib.hash_kmer(codes[i:i + kc])))
    if not km:
        return []
    seen = set()
    returned = []
    for i in range(len(km) - w1 + 1):
        best = min(km[i:i + w1], key=lambda t: t[1])
        key = (best[1], best[0] // w1)
        if key in seen:
            continue
        seen.add(key)
        returned.append(best)
    return returned


# 兼容旧名（术语统一为 anchor）
strobes = anchors


def evoseed_keys(seq: str, kc: int = 5, gmin: int = 30, gmax: int = 300,
                 delta: int = 15, use_class: bool = True, use_gap: bool = True,
                 w1: int = 6, max_pairs_per_strobe: int = 12) -> set[int]:
    """返回该序列的锚点对种子键集合（128-bit 键）。"""
    st = anchors(seq, kc, w1, use_class)
    keys = set()
    m = len(st)
    for i in range(m):
        p1, h1 = st[i]
        cnt = 0
        for j in range(i + 1, m):
            p2, h2 = st[j]
            d = p2 - p1
            if d < gmin:
                continue
            if d > gmax:
                break
            bkt = (d // delta) if use_gap else 0
            keys.add(seedlib.hash_kmer([h1, h2, bkt], salt=0xE05EED))
            cnt += 1
            if cnt >= max_pairs_per_strobe:
                break
    return keys


def build_index(sequences: dict[str, str], **kw) -> dict[int, list[str]]:
    idx: dict[int, list[str]] = {}
    for sid, s in sequences.items():
        for k in evoseed_keys(s, **kw):
            idx.setdefault(k, []).append(sid)
    return idx


def shared_score(idx_a: set[int], idx_b: set[int]) -> int:
    return len(idx_a & idx_b)


if __name__ == "__main__":
    # 自检：类空间锚点对数应远小于原始 k-mer 数，且具组合性
    demo = "MFKVLTIMGGLIALADLSAKPVGWRIRLPGGQVAEPIKVDANGVITDA" + \
           "QSSSLIFRGFERILQHQDPRDAAFLTQRICGVCPLSHGLTATNALDELYGVAEH"
    raw = evoseed_keys(demo, use_class=False)
    cls = evoseed_keys(demo, use_class=True)
    print("len(seq)=", len(demo))
    print("raw  amino-acid anchor-pair keys:", len(raw))
    print("substitution-class anchor-pair keys:", len(cls))
