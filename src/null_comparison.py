"""Two nulls for the same statistic, and why the choice decides the answer.

Ferrari et al. (Nat Commun 2019, 10:737) compared diurnal transcriptomes across
nine Archaeplastida with this statistic: the mean |phase difference| between
one-to-one orthologs of two species, optionally after shifting one species'
phases to the best offset. Their null permutes the TIME POINTS of the expression
matrix and re-runs rhythm detection, which produces near-uniform expected
phases. They concluded that diurnal programs are "remarkably similar" across the
kingdom, and themselves note that rhythmic genes pile up at dusk/night.

That pile-up is the problem. If two organisms both concentrate their rhythmic
transcripts into a few hours of the day, their orthologs will have a small mean
phase difference relative to a uniform expectation **whether or not the timing
is ortholog-specific**. A null that preserves each dataset's own phase
distribution and shuffles only which ortholog carries which phase separates the
two. This module runs both nulls on the same cyanobacterial data.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import load, rhythm, orthology as O
from analyze import OUT, phase_diff
from datasets import BY_ID

from paths import GEO, KEGG as CACHE
BRIDGE = {"MED4": "MED4.gff.gz", "HALO_05": "HAO.gff.gz", "HALO_25": "HAO.gff.gz"}
QCUT = 0.05


def prepare(did):
    """expression + metadata + a one-to-one KO map for one dataset.

    'One-to-one' here means a KO carried by exactly one gene in that genome's
    KEGG annotation and present exactly once in the dataset - the nearest
    available equivalent of the one-to-one OrthoFinder orthologs Ferrari et al.
    used, without running a full orthology pipeline.
    """
    expr, meta = load.LOADERS[did](GEO)
    d = BY_ID[did]
    g2ko = O.gene_to_ko(d.kegg, CACHE)
    if did in BRIDGE:
        bridge = O.locus_tag_bridge(GEO / BRIDGE[did])
        keys = [bridge.get(g, "") for g in expr.index]
    else:
        keys = list(expr.index)
    ko = pd.Series([g2ko.get(k, "") for k in keys], index=expr.index)
    # single-copy in the genome annotation AND unique in this dataset
    genome_counts = pd.Series(list(g2ko.values())).value_counts()
    single = set(genome_counts[genome_counts == 1].index)
    ko = ko[ko.isin(single)]
    ko = ko[~ko.duplicated(keep=False)]
    expr = expr.loc[ko.index]
    expr.index = ko.values
    return expr, meta


def phases(expr, meta, period=24.0):
    r = rhythm.cosinor(expr, meta, period)
    return r["phase_h"], r["q"]


def mean_dphase(x, y, shift=0.0):
    return float(np.mean(np.abs(phase_diff(x, y + shift))))


def best_shift(x, y, grid=np.arange(-12, 12.5, 1.0)):
    vals = [(mean_dphase(x, y, s), s) for s in grid]
    return min(vals)


def compare(a, b, n_perm=200, seed=0, use_shift=True):
    rng = np.random.default_rng(seed)
    ea, ma = prepare(a)
    eb, mb = prepare(b)
    shared = ea.index.intersection(eb.index)
    ea, eb = ea.loc[shared], eb.loc[shared]

    pa, qa = phases(ea, ma)
    pb, qb = phases(eb, mb)
    keep = (qa < QCUT) & (qb < QCUT)
    keep = keep[keep].index
    if len(keep) < 15:
        return None
    x, y = pa[keep].values, pb[keep].values
    obs, shift = best_shift(x, y) if use_shift else (mean_dphase(x, y), 0.0)

    # --- null 1 (Ferrari-style): permute time labels, re-detect rhythms, recompute
    null_time = []
    for _ in range(n_perm):
        ma2 = ma.copy(); mb2 = mb.copy()
        ma2["time"] = rng.permutation(ma2["time"].values)
        mb2["time"] = rng.permutation(mb2["time"].values)
        pa2, qa2 = phases(ea, ma2)
        pb2, qb2 = phases(eb, mb2)
        k2 = (qa2 < QCUT) & (qb2 < QCUT)
        k2 = k2[k2].index
        if len(k2) < 10:
            # with time labels destroyed few genes pass; fall back to all shared
            k2 = shared
        v, _ = best_shift(pa2[k2].values, pb2[k2].values) if use_shift else \
               (mean_dphase(pa2[k2].values, pb2[k2].values), 0)
        null_time.append(v)
    null_time = np.array(null_time)

    # --- null 2 (marginal-preserving): same phases, shuffled ortholog identity
    null_marg = []
    for _ in range(max(n_perm, 1000)):
        v, _ = best_shift(x, rng.permutation(y)) if use_shift else \
               (mean_dphase(x, rng.permutation(y)), 0)
        null_marg.append(v)
    null_marg = np.array(null_marg)

    return dict(a=a, b=b, n_orthologs=len(keep), best_shift_h=shift,
                obs_mean_dphase_h=round(obs, 2),
                null_time_mean=round(float(null_time.mean()), 2),
                p_time_null=round(float((null_time <= obs).mean()), 4),
                null_marginal_mean=round(float(null_marg.mean()), 2),
                p_marginal_null=round(float((null_marg <= obs).mean()), 4))


def main():
    pairs = [("SE7942_ito_WT", "SE7942_vijayan"),   # within-strain reference
             ("SE7942_ito_WT", "SYN6803"), ("SE7942_ito_WT", "TERY"),
             ("SE7942_ito_WT", "MED4"), ("SE7942_ito_WT", "NPUN"),
             ("SYN6803", "TERY"), ("SYN6803", "MED4"), ("SYN6803", "NPUN"),
             ("TERY", "MED4"), ("TERY", "NPUN"), ("MED4", "NPUN")]
    rows = [r for r in (compare(a, b) for a, b in pairs) if r]
    t = pd.DataFrame(rows)
    t["comparison"] = np.where((t.a == "SE7942_ito_WT") & (t.b == "SE7942_vijayan"),
                               "WITHIN-STRAIN", "cross-lineage")
    t.to_csv(OUT / "n1_two_nulls.csv", index=False)
    print(t.to_string(index=False))
    c = t[t.comparison == "cross-lineage"]
    print(f"\ncross-lineage pairs called significantly similar (p<0.05):")
    print(f"  with the time-permutation null (Ferrari-style): "
          f"{(c.p_time_null < 0.05).sum()} of {len(c)}")
    print(f"  with the marginal-preserving null:              "
          f"{(c.p_marginal_null < 0.05).sum()} of {len(c)}")


if __name__ == "__main__":
    main()
