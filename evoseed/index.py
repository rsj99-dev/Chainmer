"""EvoSeed API —— 索引层（seed_index）：倒排表构建与序列化。

索引文件包含可复现元数据：
  - format / algorithm / version
  - generate 参数（kc/gmin/gmax/delta/use_class/gap_filter/gap_bucket/w1/...）
  - 六类替换映射（class_map）
  - 哈希方案（128-bit：splitmix64 × 2）
  - 数据库信息（序列数、ID 列表 SHA256）
保存/加载后应得到相同检索结果（有单测）。
"""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from pathlib import Path

from . import generator as gen

INDEX_FORMAT = "evoseed-index/1"
ALGORITHM_VERSION = "0.2.0"
HASH_SCHEME = "splitmix64-x2/128"


class SeedIndex(dict):
    """dict 子类：键 -> postings，另带 `.meta` 元数据。"""

    def __init__(self, *args, meta=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.meta = dict(meta or {})


def _db_info(sequences: dict[str, str]) -> dict:
    h = hashlib.sha256()
    for sid in sorted(sequences):
        h.update(sid.encode("utf-8"))
        h.update(b"\n")
    return {"n_sequences": len(sequences), "ids_sha256": h.hexdigest()}


def build(sequences: dict[str, str], params: dict | None = None,
          class_map: dict | None = None, db_info: dict | None = None,
          **kw) -> SeedIndex:
    """sequences: {id: seq} -> SeedIndex（key -> [seq_id,...]），并记录元数据。

    `params` 为传给 generate() 的参数（不要与 `kw` 混用；`kw` 作为兼容别名）。
    `class_map` 建议传入 seedlib.AA_CLASS，以便索引文件记录六类映射。
    """
    p = dict(params or {})
    p.update(kw)
    effective = dict(getattr(gen, "DEFAULT_PARAMS", {}))
    effective.update(p)
    meta = {
        "format": INDEX_FORMAT,
        "algorithm": "EvoSeed",
        "version": ALGORITHM_VERSION,
        "params": effective,
        "class_map": dict(class_map or {}),
        "hash_scheme": HASH_SCHEME,
        "db": dict(db_info or _db_info(sequences)),
    }
    idx: dict[int, list[str]] = defaultdict(list)
    for sid, s in sequences.items():
        for k in gen.generate(s, **p):
            idx[k].append(sid)
    out = SeedIndex(meta=meta)
    out.update(idx)
    return out


def save(index: dict[int, list[str]], path) -> None:
    payload = {"meta": getattr(index, "meta", {}),
               "index": {str(k): v for k, v in index.items()}}
    Path(path).write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")


def load(path) -> SeedIndex:
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    if isinstance(raw, dict) and "index" in raw:
        meta = raw.get("meta", {})
        data = {int(k): v for k, v in raw["index"].items()}
    else:                     # 兼容旧版纯 dict 索引
        meta, data = {}, {int(k): v for k, v in raw.items()}
    return SeedIndex(data, meta=meta)
