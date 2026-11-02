"""EvoSeed API —— 检索层（seed_search）：对倒排索引执行候选召回。

约定：
- `search()` 只做计数召回；若要施加**负向排斥**，必须同时提供候选的键集 `db_keys`
  （因为排斥需要比较候选自身的负向键命中数），否则抛出 ValueError，绝不静默忽略。
- `max_reject`：允许命中的负向键数量上限；命中数 > max_reject 的候选被剔除。
- `search_weighted()`：使用训练折学习的 SWPN 权重做加权排序。
"""
from __future__ import annotations

from collections import Counter

from . import generator as gen
from . import scoring


def search(query_seq: str, index: dict[int, list[str]], min_shared: int = 1,
           negative_keys: set[int] | None = None, max_reject: int = 0,
           db_keys: dict[str, set[int]] | None = None,
           **params) -> list[tuple[str, int]]:
    """返回候选 [(seq_id, 共享锚点对键数)]，按分排序。

    negative_keys 非空时必须提供 db_keys（候选键集），否则无法计算排斥，直接报错。
    """
    if negative_keys and db_keys is None:
        raise ValueError("negative rejection requires db_keys (candidate key sets)")
    qk = gen.generate(query_seq, **params)
    votes: Counter = Counter()
    for k in qk:
        for sid in index.get(k, ()):
            votes[sid] += 1
    hits = [(sid, c) for sid, c in votes.items() if c >= min_shared]
    if negative_keys:
        hits = [(sid, c) for sid, c in hits
                if scoring.reject_score(db_keys[sid], negative_keys) <= max_reject]
    hits.sort(key=lambda x: -x[1])
    return hits


def search_with_keys(query_seq: str, db_keys: dict[str, set[int]],
                     min_shared: int = 1, negative_keys: set[int] | None = None,
                     max_reject: int = 0, **params):
    """候选已预生成键集时（可同时施加负向排斥），返回计数得分。"""
    return scoring.candidate_rank(query_seq, db_keys, min_shared=min_shared,
                                  negative_keys=negative_keys,
                                  max_reject=max_reject, **params)


def search_weighted(query_seq: str, db_keys: dict[str, set[int]],
                    weights: dict[int, float], min_score: float | None = None,
                    min_shared: int = 1, negative_keys: set[int] | None = None,
                    max_reject: int = 0, **params):
    """使用训练折 SWPN 权重进行加权候选排序（返回 float 得分）。"""
    return scoring.candidate_rank_weighted(
        query_seq, db_keys, weights, min_score=min_score, min_shared=min_shared,
        negative_keys=negative_keys, max_reject=max_reject, **params)
