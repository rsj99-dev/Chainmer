"""EvoSeed 单元/回归测试（纯 assert，可直接 `python tests/test_evoseed.py` 运行，
也可被 pytest 收集）。覆盖：负向排斥、SWPN 权重方向与排序、标准 syncmer、128-bit 哈希、
索引存取、边界输入。"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from evoseed import seeds as seedlib
from evoseed import core as evos_seed
from evoseed import generator as gen
from evoseed import index as sidx
from evoseed.search import (search as search_candidates, search_with_keys,
                            search_weighted_index)
from evoseed import scoring

LONG = ("MVKAVLAVLAKGGTIAAQLSPEHMPVRMRMLDNNEAAELKALGVPFGRIRLPGGQVAEPIKVDANG"
        "VITDAQSSSLIFRGFERILQHQDPRDAAFLTQRICGVCPLSHGLTATNALDELYGVAEHVPKD"
        "GYKVTYAPQISDEEVKQMLEDLKRAG")


def test_hash_deterministic():
    codes = seedlib.encode("ACDEFGH")
    assert seedlib.hash_kmer(codes) == seedlib.hash_kmer(codes)
    assert seedlib.hash_kmer(codes) != seedlib.hash_kmer(seedlib.encode("ACDEFGK"))


def test_hash_is_128bit():
    h = seedlib.hash_kmer(seedlib.encode(LONG[:20]))
    assert h.bit_length() > 64          # 128-bit 量级，碰撞可忽略
    h2 = seedlib.hash_kmer(seedlib.encode(LONG[1:21]))
    assert h != h2


def test_class_encode_lowercase():
    assert seedlib.encode("a")[0] >= 21
    cls = seedlib.to_class_string("AVLIM")
    assert set(cls) <= set("abcdef")


def test_class_mapping_disjoint():
    # 20 种标准氨基酸各自恰好属于一个类；G 只在 small_polar（b）
    seen = {}
    for aa in seedlib.AA20:
        c = seedlib.AA_CLASS[aa]
        assert c in "abcdef"
        seen.setdefault(c, set()).add(aa)
    assert seen["b"] == set("STNPG")
    assert "G" not in seen["a"]
    assert len([aa for aa in seedlib.AA20 if seedlib.AA_CLASS[aa] == "b"]) == 5


def test_clean_protein_policies():
    assert seedlib.clean_protein("ACDEFG") == "ACDEFG"
    assert seedlib.clean_protein("ACDUFG") == "ACDCFG"          # U -> C
    assert seedlib.clean_protein("ACDUFG", sec_policy="X") == "ACDXFG"
    assert seedlib.clean_protein("AC-D.F*G") == "ACDFXG"       # 空位移除, * -> X
    assert set(seedlib.clean_protein("ACBZJOXG")) <= set(seedlib.AA20 + "X")
    try:
        seedlib.clean_protein("ACD", sec_policy="Z")
    except ValueError:
        pass
    else:
        raise AssertionError("sec_policy 非法值应报错")


def test_syncmer_open_definition():
    k, s = 5, 2
    codes = seedlib.encode(LONG)
    expected = set()
    for i in range(len(codes) - k + 1):
        smers = [seedlib.hash_kmer(codes[i + j:i + j + s]) for j in range(k - s + 1)]
        mn = min(smers)
        t = smers.index(mn)                      # 平局取最左
        if t == 0 or t == k - s:
            expected.add(seedlib.hash_kmer(codes[i:i + k]))
    assert seedlib.syncmers(LONG, k, s, mode="open") == expected


def test_syncmer_closed_definition():
    k, s = 5, 2
    codes = seedlib.encode(LONG)
    expected = set()
    for i in range(len(codes) - k + 1):
        smers = [seedlib.hash_kmer(codes[i + j:i + j + s]) for j in range(k - s + 1)]
        mn = min(smers)
        if smers.index(mn) == 0:
            expected.add(seedlib.hash_kmer(codes[i:i + k]))
    assert seedlib.syncmers(LONG, k, s, mode="closed") == expected


def test_syncmer_invalid_params():
    for k, s in [(5, 5), (5, 6), (0, 1)]:
        try:
            seedlib.syncmers(LONG, k, s)
        except ValueError:
            continue
        raise AssertionError(f"k={k},s={s} 应报错")


def test_evoseed_self_match():
    s = LONG
    k = evos_seed.evoseed_keys(s)
    assert len(k & evos_seed.evoseed_keys(s)) == len(k)
    assert len(k) > 0


def test_cross_family_lower_shared():
    mvh = LONG
    frh = ("MTNKKAAAWLVGAIVGGLGMAEAKLDDKAEAKKAAKEAKKAAELKALGVPFGRIRLPGGQVAEPIKVDANG"
           "VITDAQSSSLIFRGFERILQHQDPRDAAWQISNEEVKQMLEDLKRAG")
    self_shared = len(gen.generate(mvh) & gen.generate(mvh))
    cross = len(gen.generate(mvh) & gen.generate(frh))
    assert cross < self_shared


def test_gap_ablation_changes_keys():
    with_gap = evos_seed.evoseed_keys(LONG, use_gap=True)
    no_gap = evos_seed.evoseed_keys(LONG, use_gap=False)
    assert with_gap != no_gap


def test_gap_modes_distinct():
    full = evos_seed.evoseed_keys(LONG, use_gap=True)
    filter_only = evos_seed.evoseed_keys(LONG, gap_filter=True, gap_bucket=False)
    no_gap = evos_seed.evoseed_keys(LONG, use_gap=False)
    assert len(full) > 0 and len(filter_only) > 0 and len(no_gap) > 0
    assert full != filter_only          # 去掉桶，保留距离过滤
    assert filter_only != no_gap        # 去掉距离过滤
    assert full != no_gap


def test_amino_mapping_complete():
    assert len(seedlib.AA20) == 20
    for aa in seedlib.AA20:
        assert aa in seedlib.AA_CLASS
    # 每个标准氨基酸恰好一个类；特殊字符归 f
    assert len({seedlib.AA_CLASS[a] for a in seedlib.AA20}) == 5
    for ch in "XBZJO":
        assert seedlib.AA_CLASS[ch] == "f"


def test_index_matches_brute_force():
    db = {"a": LONG, "b": LONG[:60] + "W" + LONG[61:], "c": LONG[::-1]}
    idx = sidx.build(db, class_map=seedlib.AA_CLASS)
    dbk = {k: gen.generate(v) for k, v in db.items()}
    for q in (LONG, db["b"]):
        hit_idx = dict(search_candidates(q, idx, min_shared=1))
        hit_brute = dict(search_with_keys(q, dbk, min_shared=1))
        assert hit_idx == hit_brute


def test_weighted_index_matches_brute_force():
    db = {"a": LONG, "b": LONG[:60] + "W" + LONG[61:], "c": LONG[::-1]}
    idx = sidx.build(db, class_map=seedlib.AA_CLASS)
    dbk = {k: gen.generate(v) for k, v in db.items()}
    qk = gen.generate(LONG)
    key = sorted(qk)[0]
    w = {key: 2.5}
    hit_idx = dict(search_weighted_index(LONG, idx, w, min_score=1e-9))
    hit_brute = dict(scoring.candidate_rank_weighted(LONG, dbk, w, min_shared=1,
                                                     min_score=1e-9))
    assert set(hit_idx) == set(hit_brute)
    for sid in hit_idx:
        assert abs(hit_idx[sid] - hit_brute[sid]) < 1e-9


def test_negative_rejection_filters_candidate():
    qk = gen.generate(LONG)
    keys_a = set(qk)                       # 与查询完全相同
    keys_b = set(qk)                       # 也相同，但会被负向键命中
    neg = set(list(sorted(keys_b))[:3])    # 取 3 个键作为负向集
    db = {"A": keys_a - neg, "B": keys_b}  # A 不含负向键；B 含 3 个
    hits0 = dict(scoring.candidate_rank(LONG, db, min_shared=1,
                                        negative_keys=neg, max_reject=0))
    assert "A" in hits0 and "B" not in hits0          # max_reject=0 时剔除 B
    hits2 = dict(scoring.candidate_rank(LONG, db, min_shared=1,
                                        negative_keys=neg, max_reject=2))
    assert "B" not in hits2                            # 命中 3 > 2，仍剔除
    hits3 = dict(scoring.candidate_rank(LONG, db, min_shared=1,
                                        negative_keys=neg, max_reject=3))
    assert "B" in hits3                                # 命中数 == 上限，保留（> 才剔除）


def test_search_requires_db_keys_for_rejection():
    idx = sidx.build({"A": LONG, "B": LONG})
    try:
        search_candidates(LONG, idx, negative_keys={123})
    except ValueError:
        pass
    else:
        raise AssertionError("negative_keys 但缺 db_keys 时应报错")


def test_search_rejection_with_db_keys():
    idx = sidx.build({"A": LONG, "B": LONG})
    db = {k: gen.generate(v) for k, v in {"A": LONG, "B": LONG}.items()}
    neg = set(list(sorted(db["B"]))[:3])
    db["A"] = db["A"] - neg
    hits = dict(search_candidates(LONG, idx, negative_keys=neg, max_reject=0, db_keys=db))
    assert "A" in hits and "B" not in hits


def test_contrastive_weight_direction():
    pos_inter = [{1, 2}]
    neg_inter = [{2, 3}]
    w = scoring.contrastive_weights(pos_inter, neg_inter)
    assert w[1] > 0            # 仅正样本出现
    assert abs(w[2]) < 1e-12   # 正负各一次
    assert w[3] < 0            # 仅负样本出现
    assert scoring.weighted_score({1}, {1}, w) > 0
    assert scoring.weighted_score({3}, {3}, w) < 0


def test_weighted_search_orders_positive_enriched_first():
    pos_inter = [{1, 2}]
    neg_inter = [{3}]
    w = scoring.contrastive_weights(pos_inter, neg_inter)
    db = {"pos_like": {1}, "neg_like": {3}}
    assert scoring.weighted_score({1, 9}, db["pos_like"], w) > \
           scoring.weighted_score({1, 9}, db["neg_like"], w)
    order = [sid for sid, _sc in sorted(
        [(sid, scoring.weighted_score({1, 9}, ks, w)) for sid, ks in db.items()],
        key=lambda x: -x[1])]
    assert order[0] == "pos_like"


def test_index_save_load_roundtrip():
    db = {"a": LONG, "b": LONG}
    idx = sidx.build(db, class_map=seedlib.AA_CLASS)
    assert idx.meta["algorithm"] == "EvoSeed"
    assert idx.meta["hash_scheme"].endswith("/128")
    assert idx.meta["db"]["n_sequences"] == 2
    assert idx.meta["params"]["kc"] == 5
    assert idx.meta["params"]["gap_filter"] is True
    assert idx.meta["class_map"]["G"] == "b"
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "idx.json"
        sidx.save(idx, p)
        idx2 = sidx.load(p)
    assert idx2 == idx
    assert idx2.meta == idx.meta
    r1 = search_with_keys(LONG, {k: gen.generate(v) for k, v in db.items()})
    r2 = search_with_keys(LONG, {k: gen.generate(v) for k, v in db.items()})
    assert [x[0] for x in r1] == [x[0] for x in r2]
    # 从加载后的索引检索，结果与内存索引一致
    assert dict(search_candidates(LONG, idx)) == dict(search_candidates(LONG, idx2))


def test_edge_cases():
    assert seedlib.extract_seeds("", "kmer", k=5) == set()
    assert seedlib.extract_seeds("ACD", "kmer", k=5) == set()
    assert seedlib.syncmers("ACD", 5, 2) == set()
    for s in ["ACDXBZJO*", "acdefg", "ACD--EFG"]:
        seedlib.extract_seeds(s, "kmer", k=3)      # 不应抛异常
    dup = [LONG, LONG]
    assert gen.generate(dup[0]) == gen.generate(dup[1])
    assert len(evos_seed.evoseed_keys("ACD")) == 0


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
