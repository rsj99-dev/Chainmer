"""EvoSeed API —— 打分与负向排斥层（scoring）。

两类得分：
  1) 计数得分：查询与候选共享的锚点对种子键数；
  2) 加权得分：共享键权重之和（SWPN：正负对比权重）。
负向排斥：max_reject 为**允许命中的负向键数量上限**；命中数 > max_reject 则剔除。
一次命中(共享>=1)绝不等于直系同源判定，仅用于候选排序/召回。
"""
from __future__ import annotations

from . import generator as gen


def consistency_score(keys_q: set[int], keys_c: set[int]) -> int:
    return len(keys_q & keys_c)


def reject_score(keys_c: set[int], negative_keys: set[int]) -> int:
    """候选命中了多少个负向排斥键（越大越应剔除）。"""
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
    """按**计数得分**排序的候选 [(id, score)]；命中负向键超过 max_reject 则剔除。"""
    qk = gen.generate(query_seq, **params)
    out = []
    for sid, keys in db_keys.items():
        shared = consistency_score(qk, keys)
        if shared < min_shared:
            continue
        if negative_keys:
            if reject_score(keys, negative_keys) > max_reject:
                continue
        out.append((sid, shared))
    out.sort(key=lambda x: -x[1])
    return out


def candidate_rank_weighted(query_seq: str, db_keys: dict[str, set[int]],
                            weights: dict[int, float],
                            min_score: float | None = None,
                            min_shared: int = 1,
                            negative_keys: set[int] | None = None,
                            max_reject: int = 0, **params
                            ) -> list[tuple[str, float]]:
    """按**加权得分**（SWPN）排序的候选 [(id, score)]。

    weights 必须由训练折学习；本函数不参与权重拟合。
    未在 weights 中的共享键取 default 0.0（中性）。
    """
    qk = gen.generate(query_seq, **params)
    out = []
    for sid, keys in db_keys.items():
        shared = qk & keys
        if len(shared) < min_shared:
            continue
        if negative_keys and reject_score(keys, negative_keys) > max_reject:
            continue
        sc = weighted_score(qk, keys, weights, default_w=0.0)
        if min_score is not None and sc < min_score:
            continue
        out.append((sid, sc))
    out.sort(key=lambda x: -x[1])
    return out


# --------------------------------------------------------------------------- #
# EvoSeed-SW：正交判别特异性加权（方向 C）
# 权重来自训练折中该种子键在旁系/跨亚群负样本对上的频率：
#     w(k) = log((N_neg + 1) / (freq_neg(k) + 1))
# 罕见于负样本的键（更能区分正交）获得更高权重；严格禁止用测试折学习频次。
# --------------------------------------------------------------------------- #
def specificity_weights(neg_pair_key_intersections, n_neg: int) -> dict[int, float]:
    """仅负向权重（基线）：w(k) = ln((N_neg+1)/(freq_neg(k)+1))。"""
    import math
    from collections import Counter
    freq = Counter()
    for inter in neg_pair_key_intersections:
        freq.update(inter)
    return {k: math.log((n_neg + 1) / (freq.get(k, 0) + 1)) for k in freq}


def contrastive_weights(pos_pair_key_intersections, neg_pair_key_intersections
                        ) -> dict[int, float]:
    """**SWPN 正负对比权重**：w(k) = ln((f_pos(k)+1)/(f_neg(k)+1))。

    - f_pos / f_neg 分别为训练折直系对 / 负样本对的共享键频次；
    - 仅在训练折调用；测试折不得参与。
    - 未出现在任一频次表中的键不在此字典中（加权时取 default 0.0）。
    """
    import math
    from collections import Counter
    fpos, fneg = Counter(), Counter()
    for inter in pos_pair_key_intersections:
        fpos.update(inter)
    for inter in neg_pair_key_intersections:
        fneg.update(inter)
    keys = set(fpos) | set(fneg)
    return {k: math.log((fpos.get(k, 0) + 1) / (fneg.get(k, 0) + 1)) for k in keys}


def weighted_score(keys_q: set[int], keys_c: set[int],
                   weights: dict[int, float], default_w: float = 0.0) -> float:
    """加权一致性分：共享键权重之和（SWPN 默认 0.0 为中性）。"""
    shared = keys_q & keys_c
    return sum(weights.get(k, default_w) for k in shared)
