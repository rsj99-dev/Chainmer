"""EvoSeed API —— 延伸层（seed_extend）。

对通过种子一致性的候选做精修打分。原型提供两种：
  1) alignment-free：直接使用共享锚点对键数（默认，快）；
  2) banded Smith-Waterman：Biopython PairwiseAligner(BLOSUM62) 复核（可选，精）。
延伸仅提升排序/候选质量，仍不等同于直系同源最终判定（需交系统发育层）。
"""
from __future__ import annotations

from . import generator as gen


def score_pair(query_seq: str, cand_seq: str, **params) -> int:
    """默认：共享锚点对键数。"""
    return len(gen.generate(query_seq, **params) & gen.generate(cand_seq, **params))


def score_pair_align(query_seq: str, cand_seq: str) -> float:
    """可选：全局比对位分（BLOSUM62），用于精修/验证。"""
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -11.0
    al.extend_gap_score = -1.0
    return float(al.score(query_seq, cand_seq))
