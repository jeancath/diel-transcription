"""Is the within-strain vs cross-lineage gap an artifact of double thresholding?

Brooks et al. (J Biol Rhythms 2023) showed that requiring a gene to clear a
significance threshold in each study separately manufactures apparent
non-replication, and that joint analyses recover far more consistency. The gap
reported here must therefore survive relaxing, and removing, the requirement
that an ortholog be called rhythmic in both datasets.

Effect size is reported as the percentage reduction in mean |phase difference|
relative to the marginal-preserving null, which is comparable across selection
rules and across pairs.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from analyze import OUT, phase_diff
from null_comparison import prepare, phases, best_shift, QCUT

RULES = {
    "both q<0.05":      lambda qa, qb: (qa < QCUT) & (qb < QCUT),
    "either q<0.05":    lambda qa, qb: (qa < QCUT) | (qb < QCUT),
    "all orthologs":    lambda qa, qb: pd.Series(True, index=qa.index),
    "top 200 combined": None,     # handled separately
}


def one_pair(a, b, n_perm=1000, seed=0):
    rng = np.random.default_rng(seed)
    ea, ma = prepare(a)
    eb, mb = prepare(b)
    shared = ea.index.intersection(eb.index)
    ea, eb = ea.loc[shared], eb.loc[shared]
    pa, qa = phases(ea, ma)
    pb, qb = phases(eb, mb)
    pa, qa, pb, qb = (s.loc[shared] for s in (pa, qa, pb, qb))

    out = []
    for name, rule in RULES.items():
        if rule is None:
            comb = (qa.rank() + qb.rank()).nsmallest(min(200, len(qa))).index
            keep = comb
        else:
            m = rule(qa, qb)
            keep = m[m].index
        if len(keep) < 15:
            continue
        x, y = pa[keep].values, pb[keep].values
        obs, shift = best_shift(x, y)
        null = np.array([best_shift(x, rng.permutation(y))[0] for _ in range(n_perm)])
        nm = float(null.mean())
        out.append(dict(a=a, b=b, rule=name, n=len(keep),
                        obs_mean_dphase_h=round(obs, 2), null_mean=round(nm, 2),
                        effect_pct=round((nm - obs) / nm * 100, 1),
                        p=round(float((null <= obs).mean()), 4)))
    return out


def main():
    pairs = [("SE7942_ito_WT", "SE7942_vijayan"),
             ("SE7942_ito_WT", "SYN6803"), ("SE7942_ito_WT", "TERY"),
             ("SE7942_ito_WT", "MED4"), ("SE7942_ito_WT", "NPUN"),
             ("SYN6803", "TERY"), ("SYN6803", "MED4"), ("SYN6803", "NPUN"),
             ("TERY", "MED4"), ("TERY", "NPUN"), ("MED4", "NPUN")]
    rows = []
    for a, b in pairs:
        rows += one_pair(a, b)
        print(".", end="", flush=True)
    t = pd.DataFrame(rows)
    t["comparison"] = np.where((t.a == "SE7942_ito_WT") & (t.b == "SE7942_vijayan"),
                               "WITHIN-STRAIN", "cross-lineage")
    t.to_csv(OUT / "n2_threshold_sensitivity.csv", index=False)
    print()
    piv = t.pivot_table(index="rule", columns="comparison", values="effect_pct",
                        aggfunc=["median", "max", "count"]).round(1)
    print("\n== effect size (% reduction in mean |dphase| vs marginal null)")
    print(piv.to_string())
    print("\n== per pair")
    print(t.pivot_table(index=["a", "b"], columns="rule", values="effect_pct").round(1).to_string())


if __name__ == "__main__":
    main()
