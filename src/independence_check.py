"""Does the permutation null assume something false?

The marginal-preserving null shuffles ortholog identity, which treats genes as
exchangeable. Bacterial genes are not: operons are co-transcribed, so genomically
adjacent genes tend to share a phase. If that dependence is strong, shuffling
breaks it and the null variance is too small - the test becomes anti-conservative
and small effects look more significant than they are.

This measures the dependence directly: the mean |phase difference| between
genomically adjacent gene pairs, against the same quantity for random pairs.
"""
import gzip
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from analyze import OUT, phase_diff
from null_comparison import GEO, prepare, phases, QCUT

RBH_DIR = GEO.parent / "rbh"
GENOME_FILE = {"SE7942_ito_WT": "SE7942", "SYN6803": "SYN6803",
               "TERY": "TERY", "MED4": "MED4", "NPUN": "NPUN"}


def gene_order(key):
    """original locus tag -> (contig, start), from the RefSeq GFF."""
    out = {}
    with gzip.open(RBH_DIR / f"{key}.gff.gz", "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            ol = re.search(r"old_locus_tag=([^;\n]+)", f[8])
            if ol:
                out[ol.group(1).split(",")[0]] = (f[0], int(f[3]), f[6])
    return out


def main(n_perm=400, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for did, key in GENOME_FILE.items():
        e, m = prepare(did)            # indexed by KO; need locus tags instead
        # re-derive phases on the raw dataset so genomic position is available
        import load, orthology as O
        expr, meta = load.LOADERS[did](GEO)
        if did == "MED4":
            b = O.locus_tag_bridge(GEO / "MED4.gff.gz")
            idx = [b.get(g, "") for g in expr.index]
            expr = expr.loc[[i != "" for i in idx]]
            expr.index = [i for i in idx if i != ""]
        p, q = phases(expr[~expr.index.duplicated()], meta)
        pos = gene_order(key)
        tab = pd.DataFrame({"phase": p, "q": q})
        tab = tab[[g in pos for g in tab.index]]
        tab["contig"] = [pos[g][0] for g in tab.index]
        tab["start"] = [pos[g][1] for g in tab.index]
        tab["strand"] = [pos[g][2] for g in tab.index]
        tab = tab[tab.q < QCUT].sort_values(["contig", "start"])
        if len(tab) < 100:
            continue
        # adjacent pairs on the same contig and the same strand (operon candidates)
        same = (tab.contig.values[:-1] == tab.contig.values[1:]) & \
               (tab.strand.values[:-1] == tab.strand.values[1:])
        gap = tab.start.values[1:] - tab.start.values[:-1]
        adj = same & (gap < 2000)
        d_adj = np.abs(phase_diff(tab.phase.values[:-1][adj], tab.phase.values[1:][adj]))
        ph = tab.phase.values
        d_rand = np.array([np.abs(phase_diff(ph, rng.permutation(ph))).mean()
                           for _ in range(n_perm)])
        rows.append(dict(dataset=did, n_rhythmic=len(tab), n_adjacent_pairs=int(adj.sum()),
                         mean_dphase_adjacent_h=round(float(d_adj.mean()), 2),
                         mean_dphase_random_h=round(float(d_rand.mean()), 2),
                         ratio=round(float(d_adj.mean() / d_rand.mean()), 3),
                         frac_adjacent_within_1h=round(float((d_adj <= 1).mean()), 3)))
        print(rows[-1], flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "n7_independence_check.csv", index=False)
    print("\n== Phase dependence between genomically adjacent genes")
    print(t.to_string(index=False))
    print("\nratio << 1 means adjacent genes share phase, so genes are NOT exchangeable "
          "and a gene-level permutation null understates the null variance.")


if __name__ == "__main__":
    main()
