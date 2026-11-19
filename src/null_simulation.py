"""How often does a uniform-expectation null declare two unrelated transcriptomes
"significantly similar", when both merely concentrate their peaks in part of the day?

Ferrari et al. (Nat Commun 2019) scored cross-species similarity of diurnal
programs as the mean |phase difference| between one-to-one orthologs, minimised
over a global shift, against a null built by permuting the time points of the
expression matrix and re-running rhythm detection - which expects near-uniform
phases. They also report that rhythmic genes pile up at dusk/night.

This simulation holds the truth fixed at ZERO ortholog-specific conservation:
each species' phases are drawn independently, from its own piled-up marginal
distribution, with no relationship between orthologs whatsoever. The only thing
that varies is the null used to judge the statistic. Any rejection is a false
positive by construction.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from analyze import OUT

PERIOD = 24.0
SHIFTS = np.arange(-12, 12.5, 1.0)


def stat(x, Y):
    """mean |dphase| at the best global shift, for each row of Y against x."""
    Y = np.atleast_2d(Y)
    out = np.empty(len(Y))
    for i in range(0, len(Y), 64):
        Yb = Y[i:i + 64]
        d = (x[None, None, :] - Yb[:, None, :] - SHIFTS[None, :, None]) % PERIOD
        d = np.minimum(d, PERIOD - d)
        out[i:i + 64] = d.mean(axis=2).min(axis=1)
    return out


def draw(n, frac_in_window, window_h, centre, rng):
    """`frac_in_window` of genes peak inside a window of width `window_h` centred
    on `centre`; the rest are spread uniformly over the day."""
    k = rng.random(n) < frac_in_window
    ph = np.where(k,
                  centre + (rng.random(n) - 0.5) * window_h,
                  rng.random(n) * PERIOD)
    return ph % PERIOD


def one_trial(n, frac, window_h, rng, n_null=300):
    ca, cb = rng.random() * PERIOD, rng.random() * PERIOD
    x = draw(n, frac, window_h, ca, rng)
    y = draw(n, frac, window_h, cb, rng)          # independent of x by construction
    obs = stat(x, y)[0]
    # null 1: uniform expected phases (what a time-point permutation produces)
    null_u = stat(x, rng.random((n_null, n)) * PERIOD)
    # null 2: marginal-preserving - same phase values, shuffled ortholog identity
    null_m = stat(x, np.array([rng.permutation(y) for _ in range(n_null)]))
    return float((null_u <= obs).mean()), float((null_m <= obs).mean())


def main(n_trials=120, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for n in (300, 1500):
        for frac, window_h in ((0.0, 24.0), (0.4, 6.0), (0.6, 6.0), (0.9, 6.0)):
            pu, pm = zip(*[one_trial(n, frac, window_h, rng) for _ in range(n_trials)])
            pu, pm = np.array(pu), np.array(pm)
            rows.append(dict(n_orthologs=n,
                             frac_in_6h_window=frac,
                             false_positive_rate_uniform_null=round(float((pu < .05).mean()), 3),
                             false_positive_rate_marginal_null=round(float((pm < .05).mean()), 3),
                             median_p_uniform=round(float(np.median(pu)), 4),
                             median_p_marginal=round(float(np.median(pm)), 3)))
            print(rows[-1], flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "n4_null_simulation.csv", index=False)
    print("\n== False-positive rate at alpha=0.05, truth = NO ortholog-specific conservation")
    print(t.to_string(index=False))


if __name__ == "__main__":
    main()
