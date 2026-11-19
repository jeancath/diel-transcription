"""Power across scenarios, not one scenario.

The first power analysis fixed the jitter of genuinely co-timed orthologs at
SD 2 h and used two sample sizes. A reviewer rightly asked whether the
conclusion survives tighter or looser conservation, and whether the simulated
sample sizes match the real comparisons, which range from 144 to 1690 orthologs.

The mapping from observed effect size to true conserved fraction is a property
of this simulation, not a measurement, and the spread across jitter levels is
the point: the same observed effect is compatible with quite different
combinations of how many orthologs are conserved and how tightly.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from analyze import OUT
from power_analysis import one_trial


def main(n_trials=40, seed=7):
    rng = np.random.default_rng(seed)
    rows = []
    for n in (150, 500, 1500):
        for sigma in (1.0, 2.0, 4.0):
            for c in (0.0, 0.10, 0.20, 0.40):
                res = [one_trial(n, c, sigma, 0.6, 6.0, rng, n_null=200)
                       for _ in range(n_trials)]
                p = np.array([r["p"] for r in res])
                e = np.array([r["effect_pct"] for r in res])
                rows.append(dict(n_orthologs=n, jitter_sd_h=sigma, true_conserved_frac=c,
                                 power=round(float((p < .05).mean()), 2),
                                 median_effect_pct=round(float(np.median(e)), 1)))
                print(rows[-1], flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "n12_power_grid.csv", index=False)
    print("\n== power (and median effect) by n, jitter and true conserved fraction")
    print(t.pivot_table(index=["jitter_sd_h", "true_conserved_frac"],
                        columns="n_orthologs",
                        values=["power", "median_effect_pct"]).to_string())


if __name__ == "__main__":
    main()
