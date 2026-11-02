"""EvoSeed API —— 索引层（seed_index）：倒排表构建与序列化。"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from . import generator as gen


def build(sequences: dict[str, str], **params) -> dict[int, list[str]]:
    """sequences: {id: seq} -> {seed_key: [seq_id,...]}"""
    idx: dict[int, list[str]] = defaultdict(list)
    for sid, s in sequences.items():
        for k in gen.generate(s, **params):
            idx[k].append(sid)
    return dict(idx)


def save(index: dict[int, list[str]], path) -> None:
    Path(path).write_text(json.dumps({str(k): v for k, v in index.items()}),
                          encoding="utf-8")


def load(path) -> dict[int, list[str]]:
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    return {int(k): v for k, v in raw.items()}
