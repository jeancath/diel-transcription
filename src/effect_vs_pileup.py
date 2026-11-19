"""Is the effect-size measure itself confounded by how piled-up the marginals are?

Effect size is the % reduction in mean |phase difference| relative to the
marginal-preserving null. When both lineages concentrate their transcripts into
a narrow window, that null mean is already small, so there is little room for a
percentage reduction: the measure loses dynamic range exactly where pile-up is
extreme. If effect size tracks pile-up across pairs, then the lower values among
pairs that exclude S. elongatus may reflect reduced sensitivity rather than
reduced conservation - which changes the interpretation of the
leave-one-lineage-out result.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from analyze import OUT

WITHIN = {"SE7942_ito_WT", "SE7942_vijayan"}


def main():
    eff = pd.read_csv(OUT / "n2_threshold_sensitivity.csv")
    marg = pd.read_csv(OUT / "r3_phase_marginals.csv").set_index("dataset")
    pile = marg["frac_in_busiest_6h"].to_dict()

    eff["is_within"] = eff.apply(lambda r: {r.a, r.b} == WITHIN, axis=1)
    eff["pile_a"] = eff.a.map(pile)
    eff["pile_b"] = eff.b.map(pile)
    eff["pile_mean"] = eff[["pile_a", "pile_b"]].mean(axis=1)
    eff["pile_max"] = eff[["pile_a", "pile_b"]].max(axis=1)
    eff.to_csv(OUT / "n11_effect_vs_pileup.csv", index=False)

    print("== Across cross-lineage pairs: does effect size track marginal pile-up?")
    for rule, g in eff[~eff.is_within].groupby("rule"):
        r1, p1 = stats.spearmanr(g.pile_mean, g.effect_pct)
        r2, p2 = stats.spearmanr(g.null_mean, g.effect_pct)
        print(f"  {rule:18s} n={len(g)}  rho(effect, mean pile-up)={r1:+.2f} p={p1:.3f}"
              f"   rho(effect, null mean |dphase|)={r2:+.2f} p={p2:.3f}")

    g = eff[(~eff.is_within) & (eff.rule == "both q<0.05")].copy()
    g["involves_SE"] = g.a.str.startswith("SE7942") | g.b.str.startswith("SE7942")
    print("\n== the raw numbers behind the leave-one-out result (both q<0.05)")
    print(g[["a", "b", "n", "obs_mean_dphase_h", "null_mean", "effect_pct",
             "pile_mean", "involves_SE"]]
          .sort_values("effect_pct", ascending=False).to_string(index=False))

    print("\n== absolute shift in mean |dphase| (hours), which has no dynamic-range problem")
    g["abs_gain_h"] = g.null_mean - g.obs_mean_dphase_h
    print(g.groupby("involves_SE")["abs_gain_h"].agg(["count", "median"]).round(3).to_string())
    w = eff[eff.is_within & (eff.rule == "both q<0.05")].iloc[0]
    print(f"  within-strain absolute gain: {w.null_mean - w.obs_mean_dphase_h:.2f} h")

    print("\n== excluding Prochlorococcus, the most extreme pile-up (0.91)")
    for rule, gg in eff[~eff.is_within].groupby("rule"):
        no_med = gg[~(gg.a.eq("MED4") | gg.b.eq("MED4"))]
        print(f"  {rule:18s} all {len(gg)} pairs median {gg.effect_pct.median():5.1f}%"
              f"   without MED4 ({len(no_med)} pairs) median {no_med.effect_pct.median():5.1f}%"
              f"   p<0.05: {(no_med.p < .05).sum()}/{len(no_med)}")


if __name__ == "__main__":
    main()
