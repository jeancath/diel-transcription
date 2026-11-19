"""Power surface: how much widespread timing conservation could go undetected?

Extends the grid to 60%, 80% and 100% conserved, which is what bounds the claim.
Reports Wilson 95% Monte Carlo intervals, and includes the zero-conservation row
as the size check. Outer trials and inner permutations are reported separately.

Not simulated, and therefore not excluded: conservation restricted to one
functional module rather than spread uniformly across orthologs. A program
conserved only within, say, photosynthesis would present as a small conserved
fraction here and is covered by the 10-20% rows only in the uniform-spread sense.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from simcheck import Checkpoint, cell_rng

sys.path.insert(0, str(Path(__file__).parent))
from analyze import OUT
from three_nulls_fast import stat_many, draw_marginal, mc_ci, PERIOD


def one_trial(n, c, sigma, frac, window_h, rng, n_inner):
    ca, cb = rng.random() * PERIOD, rng.random() * PERIOD
    off = rng.random() * PERIOD
    x = draw_marginal(n, frac, window_h, ca, rng)
    y = draw_marginal(n, frac, window_h, cb, rng)
    cons = rng.random(n) < c
    y = np.where(cons, (x + off + rng.normal(0, sigma, n)) % PERIOD, y)
    obs = stat_many(x, y[None, :])[0]
    null = stat_many(x, np.array([rng.permutation(y) for _ in range(n_inner)]))
    nm = float(null.mean())
    return float((null <= obs).mean()), (nm - obs) / nm * 100


def main(outer=400, inner=400, seed=23):
    ck = Checkpoint(OUT, "n14_power_surface.csv",
                    ["n_orthologs", "jitter_sd_h", "true_conserved_frac"])
    for n in (150, 500, 1500):
        for sigma in (1.0, 2.0, 4.0):
            for c in (0.0, 0.40, 0.60, 0.80, 1.00):
                if ck.done((n, sigma, c)):
                    continue
                rng = cell_rng(seed, (n, sigma, c))
                try:
                    res = [one_trial(n, c, sigma, 0.6, 6.0, rng, inner) for _ in range(outer)]
                except Exception as e:
                    ck.failed((n, sigma, c), e)
                    continue
                p = np.array([r[0] for r in res])
                e = np.array([r[1] for r in res])
                k = int((p < .05).sum())
                lo, hi = mc_ci(k, outer)
                ck.add(dict(outer_trials=outer, inner_permutations=inner,
                            n_orthologs=n, jitter_sd_h=sigma, true_conserved_frac=c,
                            power=round(k / outer, 3),
                            power_ci95=f"{lo:.3f}-{hi:.3f}",
                            median_effect_pct=round(float(np.median(e)), 1)))
    ck.finish(sort_by=["n_orthologs", "jitter_sd_h", "true_conserved_frac"],
              expected=[(n, s_, c) for n in (150, 500, 1500)
                        for s_ in (1.0, 2.0, 4.0) for c in (0.0, 0.40, 0.60, 0.80, 1.00)],
              config=dict(outer_trials=outer, inner_permutations=inner, seed=seed,
                          seeding="per-cell from (seed, n, jitter, frac_conserved)",
                          phase_concentration=0.6, window_h=6.0, alpha=0.05))


if __name__ == "__main__":
    main()
