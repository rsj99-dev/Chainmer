"""EvoSeed 命令行工具（最小可用）。"""
import argparse
import sys
from pathlib import Path

from . import __version__
from .generator import generate
from .index import build, save, load
from .search import search
from .seeds import AA_CLASS


def read_fasta(path):
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
    return seqs


def cmd_seed(args):
    seqs = read_fasta(args.input)
    total = 0
    for sid, s in seqs.items():
        n = len(generate(s))
        total += n
        print(f"{sid}\t{n}")
    print(f"# sequences={len(seqs)} total_keys={total}")
    return 0


def cmd_index(args):
    seqs = read_fasta(args.db)
    idx = build(seqs, class_map=AA_CLASS)
    save(idx, args.output)
    print(f"wrote {args.output}: keys={len(idx)} "
          f"seqs={idx.meta['db']['n_sequences']}")
    return 0


def cmd_search(args):
    db = read_fasta(args.db)
    idx = load(args.index) if args.index else build(db, class_map=AA_CLASS)
    for qid, q in read_fasta(args.query).items():
        hits = search(q, idx, min_shared=args.min_shared)
        print(f">{qid}")
        for sid, sc in hits[:args.top]:
            print(f"{sid}\t{sc}")
    return 0


def cmd_demo(args):
    db = {
        "MvhA_synthetic": "MVKAVLAVLAKGGTIAAQLSPEHMPVRMRMLDNNEAAELKALGVPFGRIRLPGGQVAEPIKVDANGVITDAQSSSLIFRGFERILQHQDPRDAAFLTQRI",
        "FrhA_synthetic": "MTNKKAAAWLVGAIVGGLGMAEAKLDDKAEAKKAAKEAKKAAELKALGVPFGRIRLPGGQVAEPIKVDANGVITDAQSSSLIFRGFERILQHQDPRDAA",
    }
    idx = build(db, class_map=AA_CLASS)
    for qid, q in db.items():
        print(qid, search(q, idx, min_shared=1)[:3])
    print("NOTE: synthetic illustrative sequences, not biological evidence.")
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(prog="evoseed", description="EvoSeed seed/index/search CLI")
    p.add_argument("--version", action="version", version=f"evoseed {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("seed", help="print seed counts for a FASTA")
    s.add_argument("-i", "--input", required=True)
    s.set_defaults(func=cmd_seed)
    s = sub.add_parser("index", help="build and save an index")
    s.add_argument("-d", "--db", required=True)
    s.add_argument("-o", "--output", required=True)
    s.set_defaults(func=cmd_index)
    s = sub.add_parser("search", help="search a query FASTA against a DB/index")
    s.add_argument("-d", "--db", required=True)
    s.add_argument("-q", "--query", required=True)
    s.add_argument("-x", "--index", default="")
    s.add_argument("--min-shared", type=int, default=1)
    s.add_argument("--top", type=int, default=5)
    s.set_defaults(func=cmd_search)
    s = sub.add_parser("demo", help="run a synthetic demo")
    s.set_defaults(func=cmd_demo)
    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
