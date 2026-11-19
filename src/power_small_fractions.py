"""The 10% and 20% rows of the power surface, with the same Monte Carlo treatment.

src/power_grid_full.py covers 0/40/60/80/100%. These two rows complete the grid
the review asked for, at identical settings so the cells are comparable.
"""
import sys
from pathlib import Path
import pandas as pd
from simcheck import Checkpoint, cell_rng
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from analyze import OUT
from power_grid_full import one_trial
from three_nulls_fast import mc_ci


def main(outer=400, inner=400, seed=29):
    ck = Checkpoint(OUT, "n17_power_small_fractions.csv",
                    ["n_orthologs", "jitter_sd_h", "true_conserved_frac"])
    for n in (150, 500, 1500):
        for sigma in (1.0, 2.0, 4.0):
            for c in (0.10, 0.20):
                if ck.done((n, sigma, c)):
                    continue
                rng = cell_rng(seed, (n, sigma, c))
                try:
                    res = [one_trial(n, c, sigma, 0.6, 6.0, rng, inner)
                           for _ in range(outer)]
                except Exception as e:
                    ck.failed((n, sigma, c), e)
                    continue
                p = np.array([r[0] for r in res]); e = np.array([r[1] for r in res])
                k = int((p < .05).sum()); lo, hi = mc_ci(k, outer)
                ck.add(dict(outer_trials=outer, inner_permutations=inner,
                            n_orthologs=n, jitter_sd_h=sigma,
                            true_conserved_frac=c, power=round(k / outer, 3),
                            power_ci95=f"{lo:.3f}-{hi:.3f}",
                            median_effect_pct=round(float(np.median(e)), 1)))
    ck.finish(sort_by=["n_orthologs", "jitter_sd_h", "true_conserved_frac"],
              expected=[(n, s_, c) for n in (150, 500, 1500)
                        for s_ in (1.0, 2.0, 4.0) for c in (0.10, 0.20)],
              config=dict(outer_trials=outer, inner_permutations=inner, seed=seed,
                          seeding="per-cell from (seed, n, jitter, frac_conserved)",
                          phase_concentration=0.6, window_h=6.0, alpha=0.05))


if __name__ == "__main__":
    main()
