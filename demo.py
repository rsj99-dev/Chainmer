"""EvoSeed minimal end-to-end demo (no external data required)."""
from evoseed import build_index, generate, search_with_keys, score_pair_align

DB = {
    "MvhA_ref": ("MVKAVLAVLAKGGTIAAQLSPEHMPVRMRMLDNNEAAELKALGVPFGRIRLPGGQVAEPI"
                 "KVDANGVITDAQSSSLIFRGFERILQHQDPRDAAFLTQRI"),
    "FrhA_ref": ("MTNKKAAAWLVGAIVGGLGMAEAKLDDKAEAKKAAKEAKKAAELKALGVPFGRIRLPGGQ"
                 "VAEPIKVDANGVITDAQSSSLIFRGFERILQHQDPRDAA"),
    "Ortholog_like": ("MVKAVLAVLAKGGTIAAQLSPEHMPVRMRMLDNNEAAELKALGVPFGRIRLPGGQV"
                      "AEPIKVDANGVITDAQSSSLIFRGFERILQHQDPRDAA"),
}

def main():
    db_keys = {sid: generate(s) for sid, s in DB.items()}
    index = build_index(DB)
    print("total seed keys:", len(index))
    query = DB["MvhA_ref"]
    hits = search_with_keys(query, db_keys, min_shared=1)
    print("query = MvhA_ref")
    for sid, shared in hits:
        al = score_pair_align(query, DB[sid]) if sid != "MvhA_ref" else None
        print(f"  {sid:14s} shared_anchor_pairs={shared:3d} "
              f"BLOSUM62={('%.0f' % al) if al else 'NA'}")
    print("A hit is a candidate only; orthology requires phylogeny.")

if __name__ == "__main__":
    main()
