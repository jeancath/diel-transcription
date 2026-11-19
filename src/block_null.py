"""A null that respects operon structure.

src/independence_check.py shows genomically adjacent same-strand genes share
phase (mean |dphase| 1.1-3.9 h against 1.6-5.9 h for random pairs; 27-59% of
adjacent pairs within 1 h). Genes are therefore NOT exchangeable, and the
gene-level permutation null used so far destroys that dependence, understating
the null variance and making small effects look more significant than they are.

Two structure-preserving nulls are implemented, both of which keep each
dataset's phase values exactly and break only the ortholog pairing:

  rotation  - order orthologs along species A's genome and circularly rotate
              species B's phase vector. Preserves the full autocorrelation of
              both vectors along that axis.
  block     - group orthologs into operon-like blocks (consecutive, same strand,
              gap < 2 kb in species A) and shuffle whole blocks, keeping genes
              within a block together.

If the small cross-lineage effects survive these, they are real. If they do not,
the honest conclusion is that they were an artifact of assuming independence.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import load, orthology as O
from analyze import OUT, phase_diff
from null_comparison import GEO, prepare, phases, QCUT
from independence_check import gene_order, GENOME_FILE

PERIOD = 24.0
SHIFTS = np.arange(-12, 12.5, 1.0)
WITHIN = ("SE7942_ito_WT", "SE7942_vijayan")


def stat(x, Y):
    Y = np.atleast_2d(Y)
    out = np.empty(len(Y))
    for i in range(0, len(Y), 64):
        Yb = Y[i:i + 64]
        d = (x[None, None, :] - Yb[:, None, :] - SHIFTS[None, :, None]) % PERIOD
        d = np.minimum(d, PERIOD - d)
        out[i:i + 64] = d.mean(axis=2).min(axis=1)
    return out


def ko_positions(did):
    """KO -> genomic position in that dataset's genome, for ordering."""
    key = GENOME_FILE.get(did if did != "SE7942_vijayan" else "SE7942_ito_WT")
    pos = gene_order(key)
    expr, _ = load.LOADERS[did](GEO)
    g2ko = O.gene_to_ko(__import__("datasets").BY_ID[did].kegg,
                        GEO.parent / "kegg")
    if did == "MED4":
        b = O.locus_tag_bridge(GEO / "MED4.gff.gz")
        tags = {g: b.get(g, "") for g in expr.index}
    else:
        tags = {g: g for g in expr.index}
    out = {}
    for g, t in tags.items():
        ko = g2ko.get(t, "")
        if ko and t in pos:
            out[ko] = pos[t]
    return out


def blocks_from_positions(order_pos, max_gap=2000):
    """consecutive, same-strand, close-together -> one block."""
    blocks, cur = [], [0]
    for i in range(1, len(order_pos)):
        (c0, s0, st0), (c1, s1, st1) = order_pos[i - 1], order_pos[i]
        if c0 == c1 and st0 == st1 and (s1 - s0) < max_gap:
            cur.append(i)
        else:
            blocks.append(cur); cur = [i]
    blocks.append(cur)
    return blocks


def compare(a, b, n_perm=1000, seed=0):
    ea, ma = prepare(a)
    eb, mb = prepare(b)
    shared = ea.index.intersection(eb.index)
    if len(shared) < 50:
        return None
    pa, qa = phases(ea.loc[shared], ma)
    pb, qb = phases(eb.loc[shared], mb)
    pos = ko_positions(a)
    kos = [k for k in shared if k in pos]
    if len(kos) < 50:
        return None
    kos = sorted(kos, key=lambda k: (pos[k][0], pos[k][1]))
    x, y = pa[kos].values, pb[kos].values
    qa_, qb_ = qa[kos].values, qb[kos].values
    keep = (qa_ < QCUT) & (qb_ < QCUT)
    if keep.sum() < 30:
        keep = np.ones(len(kos), bool)
        rule = "all orthologs"
    else:
        rule = "both q<0.05"
    xs, ys = x[keep], y[keep]
    order_pos = [pos[k] for k, m in zip(kos, keep) if m]
    n = len(xs)
    obs = stat(xs, ys)[0]
    rng = np.random.default_rng(seed)

    gene_null = stat(xs, np.array([rng.permutation(ys) for _ in range(n_perm)]))
    rot = np.array([np.roll(ys, int(k)) for k in rng.integers(1, n, n_perm)])
    rot_null = stat(xs, rot)
    blk = blocks_from_positions(order_pos)
    blk_perm = []
    for _ in range(n_perm):
        bb = [blk[i] for i in rng.permutation(len(blk))]
        blk_perm.append(ys[np.concatenate(bb)])
    blk_null = stat(xs, np.array(blk_perm))

    def pack(nl, tag):
        nm = float(nl.mean())
        return {f"effect_{tag}": round((nm - obs) / nm * 100, 1),
                f"p_{tag}": round(float((nl <= obs).mean()), 4),
                f"nullsd_{tag}": round(float(nl.std()), 3)}
    return dict(a=a, b=b, rule=rule, n=n, n_blocks=len(blk),
                obs_mean_dphase=round(obs, 2),
                **pack(gene_null, "gene"), **pack(rot_null, "rotation"),
                **pack(blk_null, "block"))


def main():
    pairs = [WITHIN, ("SE7942_ito_WT", "SYN6803"), ("SE7942_ito_WT", "TERY"),
             ("SE7942_ito_WT", "MED4"), ("SE7942_ito_WT", "NPUN")]
    rows = [r for r in (compare(a, b) for a, b in pairs) if r]
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "n8_block_null.csv", index=False)
    cols = ["a", "b", "n", "n_blocks", "effect_gene", "p_gene",
            "effect_rotation", "p_rotation", "effect_block", "p_block"]
    print(t[cols].to_string(index=False))


if __name__ == "__main__":
    main()
