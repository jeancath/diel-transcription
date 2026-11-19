"""The ten cross-lineage pairs are not ten independent observations.

They are built from five genomes, so each lineage appears in four pairs. Two
consequences are handled here:

1. Multiple testing. The pairwise p values are Benjamini-Hochberg corrected
   across the ten comparisons within each ortholog-selection rule. This is
   reported instead of the raw count of pairs below 0.05, which overstated how
   widespread the effect is.
2. Dependence on any single lineage. A leave-one-lineage-out analysis recomputes
   the summary after dropping each genome in turn, showing whether the
   conclusion rests disproportionately on one organism - in particular on
   S. elongatus, which is also the strain supplying the within-strain benchmark.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from analyze import OUT
from rhythm import bh

WITHIN = {"SE7942_ito_WT", "SE7942_vijayan"}
LINEAGE = {"SE7942_ito_WT": "S. elongatus", "SE7942_vijayan": "S. elongatus",
           "SYN6803": "Synechocystis", "TERY": "Trichodesmium",
           "MED4": "Prochlorococcus", "NPUN": "Nostoc"}


def main():
    t = pd.read_csv(OUT / "n2_threshold_sensitivity.csv")
    t["is_within"] = t.apply(lambda r: {r.a, r.b} == WITHIN, axis=1)
    cross = t[~t.is_within].copy()

    out = []
    for rule, g in cross.groupby("rule"):
        g = g.copy()
        g["q_bh"] = bh(g.p.values)
        out.append(g)
        print(f"\n== {rule}: {len(g)} cross-lineage pairs")
        print(f"   raw p<0.05: {(g.p < .05).sum()}   BH q<0.05: {(g.q_bh < .05).sum()}"
              f"   BH q<0.10: {(g.q_bh < .10).sum()}")
    cross = pd.concat(out)
    cross.to_csv(OUT / "n10_pairs_bh_corrected.csv", index=False)

    print("\n== Leave-one-lineage-out (median cross-lineage effect, pairs NOT involving X)")
    rows = []
    for rule, g in cross.groupby("rule"):
        g = g.copy()
        g["la"], g["lb"] = g.a.map(LINEAGE), g.b.map(LINEAGE)
        base = g.effect_pct.median()
        rec = dict(rule=rule, all_pairs=round(base, 1), n_all=len(g))
        for lin in sorted(set(g.la) | set(g.lb)):
            keep = g[(g.la != lin) & (g.lb != lin)]
            rec[f"drop_{lin.split()[0]}"] = (round(keep.effect_pct.median(), 1)
                                             if len(keep) else np.nan)
            rec[f"n_{lin.split()[0]}"] = len(keep)
        rows.append(rec)
    loo = pd.DataFrame(rows)
    loo.to_csv(OUT / "n10_leave_one_out.csv", index=False)
    cols = ["rule", "all_pairs"] + [c for c in loo.columns if c.startswith("drop_")]
    print(loo[cols].to_string(index=False))

    print("\n== effect size by whether the pair involves S. elongatus")
    g = cross.copy()
    g["involves_SE"] = g.a.map(LINEAGE).eq("S. elongatus") | g.b.map(LINEAGE).eq("S. elongatus")
    print(g.groupby(["rule", "involves_SE"])["effect_pct"]
           .agg(["count", "median"]).round(1).to_string())


if __name__ == "__main__":
    main()
