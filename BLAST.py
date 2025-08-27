
#!/usr/bin/env python3
"""
trinity_multi_overlap.py

Analyse overlap between two Trinity assemblies **and** their translated ORF sets.

Expected inputs (all FASTA):
    --dn_trinity   Trinity.fasta         (deÂ novo nucleotide transcripts)
    --gg_trinity   GGTrinity.fasta       (genomeâ€‘guided nucleotide transcripts)
    --dn_orf       ORF.fasta             (protein ORFs predicted from Trinity)
    --gg_orf       GGORF.fasta           (protein ORFs predicted from GGTrinity)

Pipeline per sequence class (nucleotide / protein):
1. Build BLAST+ databases for each file.
2. Perform reciprocal BLAST (BLASTN for transcripts, BLASTP for ORFs).
3. Filter hits by percent identity (â€‘â€‘identity) and bidirectional coverage (â€‘â€‘coverage).
4. Identify reciprocal bestâ€‘hit (RBH) pairs â†’ **intersection**.
5. Report counts of **unique** and **shared** elements; export ID lists; draw optional Venn diagrams.

Outputs (<outdir>/):
    transcripts/summary.tsv
    transcripts/intersection.txt   (Qid\tSid RBH pairs)
    transcripts/dn_unique.txt
    transcripts/gg_unique.txt

    orfs/summary.tsv
    orfs/intersection.txt          (Qid\tSid RBH pairs)
    orfs/dn_unique.txt
    orfs/gg_unique.txt

    venn_transcripts.png            (if matplotlibâ€‘venn installed)
    venn_orfs.png                   (if matplotlibâ€‘venn installed)

Dependencies:
    BLAST+ (makeblastdb, blastn, blastp)
    Python â‰¥3.8, Biopython
    matplotlib & matplotlibâ€‘venn  (optional, for Venn plots)

Example:
    python trinity_multi_overlap.py \
        --dn_trinity Trinity.fasta \
        --gg_trinity GGTrinity.fasta \
        --dn_orf ORF.fasta \
        --gg_orf GGORF.fasta \
        --identity 95 \
        --coverage 0.8 \
        --threads 8 \
        --outdir overlap4
"""

import argparse
import subprocess
import sys
from pathlib import Path
from collections import defaultdict
import csv
import textwrap

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("Biopython not found. Install with: pip install biopython")

try:
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2
    HAS_VENN = True
except ImportError:
    HAS_VENN = False

def run(cmd: str):
    """Run shell command, exit on failure."""
    print(f"[cmd] {cmd}", file=sys.stderr)
    res = subprocess.run(cmd, shell=True)
    if res.returncode != 0:
        sys.exit(f"Command failed (exit {res.returncode}): {cmd}")


def make_blast_db(fasta: Path, db_prefix: Path, dbtype: str):
    if not (db_prefix.with_suffix(".nsq") if dbtype == "nucl" else db_prefix.with_suffix(".psq")).exists():
        run(f"makeblastdb -in {fasta} -dbtype {dbtype} -out {db_prefix}")
    else:
        print(f"[info] BLAST DB for {fasta.name} exists; skipping", file=sys.stderr)


def run_blast(query: Path, db_prefix: Path, out_tsv: Path, threads: int, seqtype: str):
    fmt = "6 qseqid sseqid pident length qlen slen"
    prog = "blastn" if seqtype == "nucl" else "blastp"
    extra = "-word_size 11 -dust no" if seqtype == "nucl" else "-seg no"
    cmd = textwrap.dedent(f"""\
    {prog} -query {query} -db {db_prefix} \
          -out {out_tsv} -outfmt '{fmt}' \
          -evalue 1e-20 {extra} \
          -max_target_seqs 10 -num_threads {threads}
    """)
    run(cmd)


def parse_hits(tsv: Path, pid_thr: float, cov_thr: float):
    """Return (qualifying pair set, bestâ€‘hit dict)."""
    pairs = set()
    best = defaultdict(lambda: (0.0, ""))  # score, subj
    with tsv.open() as fh:
        for line in fh:
            q, s, pident, aln_len, qlen, slen = line.rstrip().split("\t")
            pident = float(pident)
            aln_len = int(aln_len)
            qlen = int(qlen)
            slen = int(slen)
            if pident < pid_thr:
                continue
            qcov = aln_len / qlen
            scov = aln_len / slen
            if qcov >= cov_thr and scov >= cov_thr:
                pairs.add((q, s))
                score = pident * min(qcov, scov)
                if score > best[q][0]:
                    best[q] = (score, s)
    return pairs, {q: s for q, (score, s) in best.items()}


def reciprocal_best(dn_best: dict, gg_best: dict):
    return {(q, s) for q, s in dn_best.items() if gg_best.get(s) == q}


def fasta_ids(fasta: Path):
    return {rec.id.split(" ")[0] for rec in SeqIO.parse(fasta, "fasta")}


def write_ids(path: Path, items):
    with path.open("w") as fh:
        for x in sorted(items):
            print(x, file=fh)


def venn_plot(a_only: int, b_only: int, both: int, labels: tuple[str, str], out_png: Path):
    if not HAS_VENN:
        print("[warn] matplotlibâ€‘venn missing; skip plotting", file=sys.stderr)
        return
    plt.figure(figsize=(4, 4))
    venn2(subsets=(a_only, b_only, both), set_labels=labels)
    plt.title("Overlap ({} vs {})".format(*labels))
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def analyse_pair(dn_fa: Path, gg_fa: Path, outdir: Path, seqtype: str, pid_thr: float, cov_thr: float, threads: int):
    outdir.mkdir(parents=True, exist_ok=True)
    dn_db = outdir / "dn_db"
    gg_db = outdir / "gg_db"

    dbtype = "nucl" if seqtype == "nucl" else "prot"
    make_blast_db(dn_fa, dn_db, dbtype)
    make_blast_db(gg_fa, gg_db, dbtype)

    dn_vs_gg = outdir / "dn_vs_gg.tsv"
    gg_vs_dn = outdir / "gg_vs_dn.tsv"

    run_blast(dn_fa, gg_db, dn_vs_gg, threads, seqtype)
    run_blast(gg_fa, dn_db, gg_vs_dn, threads, seqtype)

    _, dn_best = parse_hits(dn_vs_gg, pid_thr, cov_thr)
    _, gg_best = parse_hits(gg_vs_dn, pid_thr, cov_thr)

    rbhs = reciprocal_best(dn_best, gg_best)
    dn_ids = fasta_ids(dn_fa)
    gg_ids = fasta_ids(gg_fa)

    dn_inter = {q for q, _ in rbhs}
    gg_inter = {s for _, s in rbhs}

    dn_unique = dn_ids - dn_inter
    gg_unique = gg_ids - gg_inter

    # write lists
    write_ids(outdir / "intersection.txt", [f"{q}\t{s}" for q, s in rbhs])
    write_ids(outdir / "dn_unique.txt", dn_unique)
    write_ids(outdir / "gg_unique.txt", gg_unique)

    with (outdir / "summary.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Category", "Count"])
        w.writerow(["de_novo_only", len(dn_unique)])
        w.writerow(["genome_guided_only", len(gg_unique)])
        w.writerow(["intersection", len(rbhs)])

    venn_plot(len(dn_unique), len(gg_unique), len(rbhs), labels=("deÂ novo", "genomeâ€‘guided"),
              out_png=outdir.parent / f"venn_{outdir.name}.png")

    print(f"[done] {seqtype} results in {outdir}", file=sys.stderr)


    p = argparse.ArgumentParser(description="Analyse overlap between Trinity assemblies and their ORFs (nucleotide & protein levels).")
    p.add_argument("--dn_trinity", required=True, type=Path)
    p.add_argument("--gg_trinity", required=True, type=Path)
    p.add_argument("--dn_orf", required=True, type=Path)
    p.add_argument("--gg_orf", required=True, type=Path)
    p.add_argument("--identity", type=float, default=95, help="percent identity threshold (default 95)")
    p.add_argument("--coverage", type=float, default=0.8, help="bidirectional coverage threshold (default 0.8)")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--outdir", default="overlap4", type=Path)
    args = p.parse_args()

    pid_thr = args.identity
    cov_thr = args.coverage

    # transcripts (nucleotide)
    analyse_pair(args.dn_trinity, args.gg_trinity, args.outdir / "transcripts", "nucl",
                 pid_thr, cov_thr, args.threads)

    # ORFs (protein)
    analyse_pair(args.dn_orf, args.gg_orf, args.outdir / "orfs", "prot",
                 pid_thr, cov_thr, args.threads)

    print("[all done] Results summarised in", args.outdir, file=sys.stderr)


if __name__ == "__main__":
    main()

