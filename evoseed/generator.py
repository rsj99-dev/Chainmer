"""EvoSeed API —— 种子生成层（seed_generator）。

把 algorithm/evos_seed.py 的核心封装为清晰的生成接口，对应任务书第六阶段目录结构中的
seed_generator.py。
"""
from __future__ import annotations

from . import core as evos_seed

DEFAULT_PARAMS = dict(kc=5, gmin=30, gmax=300, delta=15,
                      use_class=True, use_gap=True, w1=6)


def generate(seed: str, **params) -> set[int]:
    """由一条蛋白序列生成 EvoSeed 锚点对种子键集合。"""
    kw = dict(DEFAULT_PARAMS)
    kw.update(params)
    return evos_seed.evoseed_keys(seed, **kw)


def generate_from_fasta(path, **params) -> dict[str, set[int]]:
    """读取 FASTA 并生成每条序列的种子键集合。"""
    from pathlib import Path
    seqs, name, buf = {}, None, []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        elif line.strip():
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return {sid: generate(s, **params) for sid, s in seqs.items()}
