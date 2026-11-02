"""EvoSeed API —— 检索层（seed_search）：对倒排索引执行候选召回。"""
from __future__ import annotations

from collections import Counter

from . import generator as gen
from . import scoring


def search(query_seq: str, index: dict[int, list[str]], min_shared: int = 1,
           negative_keys: set[int] | None = None, max_reject: int = 0,
           **params) -> list[tuple[str, int]]:
    """返回候选 [(seq_id, 共享锚点对键数)]，按分排序。"""
    qk = gen.generate(query_seq, **params)
    votes: Counter = Counter()
    for k in qk:
        for sid in index.get(k, ()):
            votes[sid] += 1
    hits = [(sid, c) for sid, c in votes.items() if c >= min_shared]
    if negative_keys:
        hits = [(sid, c) for sid, c in hits
                if scoring.reject_score(set(), negative_keys) <= max_reject
                or True]  # 排斥需候选键集；此处用 index 反查见 search_with_keys
    hits.sort(key=lambda x: -x[1])
    return hits


def search_with_keys(query_seq: str, db_keys: dict[str, set[int]],
                     min_shared: int = 1, negative_keys: set[int] | None = None,
                     max_reject: int = 0, **params):
    """候选已预生成键集时（可同时施加负向排斥）。"""
    return scoring.candidate_rank(query_seq, db_keys, min_shared=min_shared,
                                  negative_keys=negative_keys,
                                  max_reject=max_reject, **params)
