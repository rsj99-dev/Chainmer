"""EvoSeed API —— 打分与负向排斥层（scoring）。

一致性分 = 查询与候选共享的锚点对种子键数；负向排斥 = 命中学习到的"旁系特异键集"则降权。
一次命中(共享>=1)绝不等于直系同源判定，仅用于候选排序/召回。
"""
from __future__ import annotations

from . import generator as gen


def consistency_score(keys_q: set[int], keys_c: set[int]) -> int:
    return len(keys_q & keys_c)


def reject_score(keys_c: set[int], negative_keys: set[int]) -> int:
    """候选命中了多少个负向排斥键（越大越应丢弃）。"""
    return len(keys_c & negative_keys)


def build_negative_set(paralog_pairs_keys) -> set[int]:
    """从旁系/近缘配对中，收集在旁系间频繁共享的锚点对键作为排斥集。

    paralog_pairs_keys: iterable of (keysA & keysB) 交集集合。
    """
    from collections import Counter
    freq = Counter()
    for inter in paralog_pairs_keys:
        freq.update(inter)
    # 至少出现 2 次的键进入排斥集（阈值可调）
    return {k for k, c in freq.items() if c >= 2}


def candidate_rank(query_seq: str, db_keys: dict[str, set[int]],
                   min_shared: int = 1, negative_keys: set[int] | None = None,
                   max_reject: int = 0, **params) -> list[tuple[str, int]]:
    """返回按一致性分排序的候选 [(id, score)]，剔除负向排斥超标者。"""
    qk = gen.generate(query_seq, **params)
    out = []
    for sid, keys in db_keys.items():
        shared = consistency_score(qk, keys)
        if shared < min_shared:
            continue
        if negative_keys:
            rej = reject_score(keys, negative_keys)
            if rej > max_reject:
                continue
        out.append((sid, shared))
    out.sort(key=lambda x: -x[1])
    return out


# --------------------------------------------------------------------------- #
# EvoSeed-SW：正交判别特异性加权（方向 C）
# 权重来自训练折中该种子键在旁系/跨亚群负样本对上的频率：
#     w(k) = log((N_neg + 1) / (freq_neg(k) + 1))
# 罕见于负样本的键（更能区分正交）获得更高权重；严格禁止用测试折学习频次。
# --------------------------------------------------------------------------- #
def specificity_weights(neg_pair_key_intersections, n_neg: int) -> dict[int, float]:
    """由训练折负样本对的"共享键集合"列表学习权重。"""
    import math
    from collections import Counter
    freq = Counter()
    for inter in neg_pair_key_intersections:
        freq.update(inter)
    return {k: math.log((n_neg + 1) / (freq.get(k, 0) + 1)) for k in freq}


def weighted_score(keys_q: set[int], keys_c: set[int],
                   weights: dict[int, float], default_w: float = 1.0) -> float:
    """特异性加权一致性分：共享键权重之和。"""
    shared = keys_q & keys_c
    return sum(weights.get(k, default_w) for k in shared)
