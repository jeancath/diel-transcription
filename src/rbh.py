"""Reciprocal best hits between proteomes - one-to-one orthologs without KEGG.

KEGG Orthology is a curated, many-to-one annotation: convenient, but it is not
an orthology assignment between two particular genomes, and KOs absent from the
annotation are silently dropped. Ferrari et al. used one-to-one OrthoFinder
orthologs; this module computes the equivalent directly with Smith-Waterman
reciprocal best hits (pyswrd), so the cross-lineage result does not rest on
KEGG's coverage.

Dataset gene ids are the ORIGINAL locus tags (Synpcc7942_0001, sll0001,
Npun_F0001, Tery_0001, PMM0001); current RefSeq proteomes carry the modern RS
tags, so the GFF old_locus_tag attribute bridges the two.
"""
import gzip
import pickle
import re
from pathlib import Path

from paths import RBH as RBH_DIR
GENOMES = ["SE7942", "SYN6803", "NPUN", "TERY", "MED4"]


def _faa(path):
    out, name, buf = {}, None, []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def _old_tag_to_protein(key):
    """original locus tag -> protein accession, via GFF old_locus_tag + feature table."""
    lt2old = {}
    with gzip.open(RBH_DIR / f"{key}.gff.gz", "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            lt = re.search(r"locus_tag=([^;\n]+)", f[8])
            ol = re.search(r"old_locus_tag=([^;\n]+)", f[8])
            if lt and ol:
                lt2old[lt.group(1)] = ol.group(1).split(",")[0]
    out = {}
    with gzip.open(RBH_DIR / f"{key}.ft.gz", "rt", errors="replace") as fh:
        head = next(fh).rstrip("\n").lstrip("# ").split("\t")
        ip, il, ic = head.index("product_accession"), head.index("locus_tag"), head.index("# feature") if "# feature" in head else 0
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(ip, il) or f[0] != "CDS" or not f[ip]:
                continue
            old = lt2old.get(f[il])
            if old:
                out[old] = f[ip]
    return out


def proteome(key):
    """original locus tag -> protein sequence."""
    seqs = _faa(RBH_DIR / f"{key}.faa.gz")
    bridge = _old_tag_to_protein(key)
    return {old: seqs[acc] for old, acc in bridge.items() if acc in seqs}


# Everything that determines an RBH call, recorded so the result can be
# reproduced or contested. pyswrd exposes Smith-Waterman via SWORD; the
# scoring matrix and gap penalties are its defaults, which are BLOSUM62 with
# gap open 10 and extend 2. Ties on score are resolved by keeping the hit
# already stored, i.e. the first encountered in pyswrd's output order - a
# strict ">" comparison below, not ">=". Ties are rare but the rule is stated
# because it is otherwise invisible.
RBH_PARAMS = dict(
    tool="pyswrd (SWORD heuristic + Smith-Waterman)",
    max_evalue=1e-5,
    max_alignments=1,
    max_candidates=30000,
    scoring_matrix="BLOSUM62 (pyswrd default)",
    gap_open=10,
    gap_extend=2,
    tie_rule="first hit encountered wins; strictly greater score replaces",
    reciprocity="forward best hit of a in b must have a as its best hit in a",
    sequence_source="NCBI RefSeq protein FASTA per assembly",
    locus_bridge="RefSeq GFF locus_tag, with old_locus_tag as fallback",
)


def _best_hits(q, t, evalue=1e-5, max_candidates=30000):
    """Best hit per query.

    NOTE: pyswrd's max_candidates behaves as a cap across the search, not per
    query - lowering it to 1000 cut Synechococcus x Nostoc from ~1500 expected
    reciprocal hits to 157. It is left at the default.
    """
    import pyswrd
    qk, qs = list(q.keys()), list(q.values())
    tk, ts = list(t.keys()), list(t.values())
    best = {}
    for hit in pyswrd.search(qs, ts, max_alignments=1, max_evalue=evalue,
                             max_candidates=max_candidates, threads=0):
        qi, ti, sc = hit.query_index, hit.target_index, hit.score
        if qk[qi] not in best or sc > best[qk[qi]][1]:
            best[qk[qi]] = (tk[ti], sc)
    return {k: v[0] for k, v in best.items()}


def rbh(a, b, cache=True):
    """Reciprocal best hits between two genomes: {locus_tag_a: locus_tag_b}."""
    f = RBH_DIR / f"rbh_{a}__{b}.pkl"
    if cache and f.exists():
        return pickle.loads(f.read_bytes())
    pa, pb = proteome(a), proteome(b)
    fwd = _best_hits(pa, pb)
    rev = _best_hits(pb, pa)
    pairs = {x: y for x, y in fwd.items() if rev.get(y) == x}
    if cache:
        f.write_bytes(pickle.dumps(pairs))
    return pairs


# The four pairs run by default span the divergence range while avoiding the
# 6900-protein Nostoc proteome, whose all-vs-all searches dominate the runtime.
SUBSET = [("SE7942", "SYN6803"), ("SE7942", "TERY"),
          ("SE7942", "MED4"), ("TERY", "MED4")]

if __name__ == "__main__":
    import itertools, sys
    pairs = SUBSET if "--subset" in sys.argv else list(itertools.combinations(GENOMES, 2))
    for a, b in pairs:
        pr = rbh(a, b)
        print(f"{a:8s} x {b:8s} {len(pr):5d} reciprocal best hits", flush=True)
