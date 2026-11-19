"""Three nulls at publication-grade replication, with Monte Carlo intervals.

src/three_nulls.py established the comparison at 50 outer trials per cell. That
is too few to quote a false-positive rate to two decimals: at p=0.88 the Monte
Carlo standard error is ~0.046, so the apparent rise from 0.76 to 0.88 with
sample size was not separable from noise.

This is the same experiment with the per-permutation cosinor fit rewritten in
pure numpy - no DataFrame construction inside the inner loop - which makes 1000
outer trials affordable. Two replication counts are reported separately, as they
answer different questions:

  OUTER trials   independent simulated experiments; sets the precision of the
                 false-positive RATE
  INNER permutations  permutations within one experiment; sets the resolution of
                 that experiment's p value

Truth is fixed at zero ortholog-specific conservation, so every rejection is a
false positive by construction.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from simcheck import Checkpoint, cell_rng

sys.path.insert(0, str(Path(__file__).parent))
from analyze import OUT

PERIOD = 24.0
SHIFTS = np.arange(-12, 12.5, 1.0)
COS_S, SIN_S = np.cos(2 * np.pi * SHIFTS / PERIOD), np.sin(2 * np.pi * SHIFTS / PERIOD)


def stat_many(x, Y, block=128):
    """mean |dphase| at the best global shift, for each row of Y against x."""
    Y = np.atleast_2d(Y)
    out = np.empty(len(Y))
    for i in range(0, len(Y), block):
        Yb = Y[i:i + block]
        d = (x[None, None, :] - Yb[:, None, :] - SHIFTS[None, :, None]) % PERIOD
        np.minimum(d, PERIOD - d, out=d)
        out[i:i + block] = d.mean(axis=2).min(axis=1)
    return out


def design(times):
    t = np.asarray(times, float)
    return np.column_stack([np.ones_like(t), np.cos(2 * np.pi * t / PERIOD),
                            np.sin(2 * np.pi * t / PERIOD)])


def phases_np(Y, times):
    """cosinor peak phase for every row of Y, in pure numpy."""
    X = design(times)
    coef = np.linalg.pinv(X) @ Y.T           # 3 x genes
    return (np.arctan2(coef[2], coef[1]) * PERIOD / (2 * np.pi)) % PERIOD


def draw_marginal(n, frac, window_h, centre, rng):
    k = rng.random(n) < frac
    return np.where(k, centre + (rng.random(n) - 0.5) * window_h,
                    rng.random(n) * PERIOD) % PERIOD


def one_trial(n, frac, window_h, times, rng, n_inner, amp=1.0, noise=0.5):
    ca, cb = rng.random() * PERIOD, rng.random() * PERIOD
    pa = draw_marginal(n, frac, window_h, ca, rng)
    pb = draw_marginal(n, frac, window_h, cb, rng)      # independent of pa
    t = np.asarray(times, float)
    Ya = amp * np.cos(2 * np.pi * (t[None, :] - pa[:, None]) / PERIOD) + \
        rng.normal(0, noise, (n, len(t)))
    Yb = amp * np.cos(2 * np.pi * (t[None, :] - pb[:, None]) / PERIOD) + \
        rng.normal(0, noise, (n, len(t)))
    fa, fb = phases_np(Ya, t), phases_np(Yb, t)
    obs = stat_many(fa, fb[None, :])[0]

    nu = stat_many(fa, rng.random((n_inner, n)) * PERIOD)

    nf = np.empty(n_inner)                               # time-point permutation
    for i in range(n_inner):
        ta, tb = rng.permutation(t), rng.permutation(t)
        nf[i] = stat_many(phases_np(Ya, ta), phases_np(Yb, tb)[None, :])[0]

    nm = stat_many(fa, np.array([rng.permutation(fb) for _ in range(n_inner)]))
    return ((nu <= obs).mean(), (nf <= obs).mean(), (nm <= obs).mean(),
            nu.mean(), nf.mean(), nm.mean())


def mc_ci(k, n):
    """Wilson interval for a Monte Carlo rate."""
    if n == 0:
        return (np.nan, np.nan)
    z, p = 1.96, k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0, c - h), min(1, c + h))


def main(outer=1000, inner=200, seed=11):
    times = np.arange(0, 48, 2.0)
    ck = Checkpoint(OUT, "n13_three_nulls_1000.csv", ["n_orthologs", "frac_in_6h_window"])
    for n in (300, 1000):
        for frac in (0.0, 0.4, 0.6, 0.9):
            if ck.done((n, frac)):
                continue
            rng = cell_rng(seed, (n, frac))
            try:
                res = [one_trial(n, frac, 6.0, times, rng, inner) for _ in range(outer)]
            except Exception as e:
                ck.failed((n, frac), e)
                continue
            pu, pf, pm, mu, mf, mm = map(np.array, zip(*res))
            rec = dict(outer_trials=outer, inner_permutations=inner,
                       n_orthologs=n, frac_in_6h_window=frac)
            for tag, p in (("uniform", pu), ("timepoint", pf), ("marginal", pm)):
                k = int((p < .05).sum())
                lo, hi = mc_ci(k, outer)
                rec[f"FPR_{tag}"] = round(k / outer, 3)
                rec[f"FPR_{tag}_ci95"] = f"{lo:.3f}-{hi:.3f}"
            rec.update(null_mean_uniform_h=round(float(mu.mean()), 2),
                       null_mean_timepoint_h=round(float(mf.mean()), 2),
                       null_mean_marginal_h=round(float(mm.mean()), 2))
            ck.add(rec)
    ck.finish(sort_by=["n_orthologs", "frac_in_6h_window"],
              expected=[(n, f) for n in (300, 1000) for f in (0.0, 0.4, 0.6, 0.9)],
              config=dict(outer_trials=outer, inner_permutations=inner, seed=seed,
                          seeding="per-cell from (seed, n, frac)",
                          times="0-46 h, 2 h steps", window_h=6.0,
                          amp=1.0, noise=0.5, alpha=0.05))


if __name__ == "__main__":
    main()
