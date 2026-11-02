"""EvoSeed: substitution-class ordered anchor-pair seeds.

Minimal public API. Pure Python (standard library only).
"""
from .core import evoseed_keys
from .generator import generate, generate_from_fasta, DEFAULT_PARAMS
from .index import build as build_index, save as save_index, load as load_index
from .search import search, search_with_keys, search_weighted, search_weighted_index
from .extend import score_pair, score_pair_align
from .scoring import (consistency_score, weighted_score, specificity_weights,
                      contrastive_weights, build_negative_set, candidate_rank,
                      candidate_rank_weighted, reject_score)

__version__ = "0.2.0"

__all__ = [
    "evoseed_keys", "generate", "generate_from_fasta", "DEFAULT_PARAMS",
    "build_index", "save_index", "load_index", "search", "search_with_keys",
    "search_weighted", "search_weighted_index", "score_pair", "score_pair_align",
    "consistency_score", "weighted_score", "specificity_weights",
    "contrastive_weights", "build_negative_set", "candidate_rank",
    "candidate_rank_weighted", "reject_score", "__version__",
]
