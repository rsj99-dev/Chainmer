"""EvoSeed 单元/回归测试（纯 assert，可直接 `python tests/test_evoseed.py` 运行，
也可被 pytest 收集）。"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from evoseed import seeds as seedlib
from evoseed import core as evos_seed
from evoseed import generator as gen
from evoseed import index as sidx
from evoseed.search import search_with_keys


def test_hash_deterministic():
    codes = seedlib.encode("ACDEFGH")
    assert seedlib.hash_kmer(codes) == seedlib.hash_kmer(codes)
    assert seedlib.hash_kmer(codes) != seedlib.hash_kmer(seedlib.encode("ACDEFGK"))


def test_class_encode_lowercase():
    # 等价类小写字母应有独立编码，不被当作 X
    assert seedlib.encode("a")[0] >= 21
    cls = seedlib.to_class_string("AVLIM")
    assert set(cls) <= set("abcdef")


def test_evoseed_self_match():
    s = "MVKAVLAVLAKGGTIAAQLSPEHMPVRMRMLDNNEAAELKALGVPFGRIRLPGGQVAEPI"
    k = evos_seed.evoseed_keys(s)
    assert len(k & evos_seed.evoseed_keys(s)) == len(k)   # 自比全共享
    assert len(k) > 0


def test_cross_family_lower_shared():
    mvh = ("MVKAVLAVLAKGGTIAAQLSPEHMPVRMRMLDNNEAAELKALGVPFGRIRLPGGQVAEPIKVDANG"
           "VITDAQSSSLIFRGFERILQHQDPRDAAFLTQRI")
    frh = ("MTNKKAAAWLVGAIVGGLGMAEAKLDDKAEAKKAAKEAKKAAELKALGVPFGRIRLPGGQVAEPIKVDANG"
           "VITDAQSSSLIFRGFERILQHQDPRDAA")
    self_shared = len(gen.generate(mvh) & gen.generate(mvh))
    cross = len(gen.generate(mvh) & gen.generate(frh))
    assert cross < self_shared


def test_gap_ablation_changes_keys():
    s = "MVKAVLAVLAKGGTIAAQLSPEHMPVRMRMLDNNEAAELKALGVPFGRIRLPGGQVAEPIKVDANGVITDA"
    with_gap = evos_seed.evoseed_keys(s, use_gap=True)
    no_gap = evos_seed.evoseed_keys(s, use_gap=False)
    assert with_gap != no_gap      # 间距分桶确为键的一部分


def test_index_roundtrip():
    long = ("MVKAVLAVLAKGGTIAAQLSPEHMPVRMRMLDNNEAAELKALGVPFGRIRLPGGQVAEPIKVDANG"
            "VITDAQSSSLIFRGFERILQHQDPRDAAFLTQRICGVCPLSHGLTATNALDELYGVAEHVPKD")
    db = {"a": long, "b": long, "c": long[:40] + "X" * 60}   # a==b 完全相同
    db_keys = {k: gen.generate(v) for k, v in db.items()}
    assert db_keys["a"], "长序列应产生非空锚点对键集"
    hits = dict(search_with_keys(long, db_keys, min_shared=1))
    assert hits["a"] == hits["b"] >= hits["c"]   # 相同序列共享最多
    assert hits["a"] > 1


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
        passed += 1
    print(f"\n{passed}/{len(tests)} tests passed")


if __name__ == "__main__":
    main()
