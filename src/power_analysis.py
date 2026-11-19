"""What size of true conservation could this design actually detect?

The reported cross-lineage effect is a 4-6% reduction in mean |phase difference|
against the marginal-preserving null. That number is meaningless on its own: a
reviewer will ask what biological truth it corresponds to, and whether the design
could have detected a moderate effect had one existed.

So here the truth is KNOWN. A fraction `c` of orthologs genuinely keep their
timing between the two species - species B's phase equals species A's plus a
fixed offset plus jitter of SD `sigma` - while the remaining 1-c are drawn
independently from the same piled-up marginal distribution. Sweeping c gives
(a) the power curve, and (b) the mapping from observed effect size back to the
fraction of genuinely co-timed orthologs.
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
    Y = np.atleast_2d(Y)
    out = np.empty(len(Y))
    for i in range(0, len(Y), 64):
        Yb = Y[i:i + 64]
        d = (x[None, None, :] - Yb[:, None, :] - SHIFTS[None, :, None]) % PERIOD
        d = np.minimum(d, PERIOD - d)
        out[i:i + 64] = d.mean(axis=2).min(axis=1)
    return out


def draw_marginal(n, frac, window_h, centre, rng):
    k = rng.random(n) < frac
    return np.where(k, centre + (rng.random(n) - 0.5) * window_h,
                    rng.random(n) * PERIOD) % PERIOD


def one_trial(n, c, sigma, frac, window_h, rng, n_null=300):
    """c = fraction of orthologs whose timing is genuinely conserved."""
    ca, cb = rng.random() * PERIOD, rng.random() * PERIOD
    offset = rng.random() * PERIOD
    x = draw_marginal(n, frac, window_h, ca, rng)
    y = draw_marginal(n, frac, window_h, cb, rng)          # the unconserved part
    conserved = rng.random(n) < c
    y = np.where(conserved, (x + offset + rng.normal(0, sigma, n)) % PERIOD, y)
    obs = stat(x, y)[0]
    null = stat(x, np.array([rng.permutation(y) for _ in range(n_null)]))
    nm = float(null.mean())
    return dict(p=float((null <= obs).mean()),
                effect_pct=(nm - obs) / nm * 100)


def main(n_trials=60, seed=1):
    rng = np.random.default_rng(seed)
    rows = []
    for n in (300, 1500):
        for c in (0.0, 0.05, 0.10, 0.20, 0.40, 0.70, 1.0):
            res = [one_trial(n, c, 2.0, 0.6, 6.0, rng) for _ in range(n_trials)]
            p = np.array([r["p"] for r in res])
            e = np.array([r["effect_pct"] for r in res])
            rows.append(dict(n_orthologs=n, true_conserved_frac=c,
                             power_at_0_05=round(float((p < .05).mean()), 3),
                             median_effect_pct=round(float(np.median(e)), 1),
                             effect_iqr=f"{np.percentile(e,25):.1f}-{np.percentile(e,75):.1f}"))
            print(rows[-1], flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "n6_power_analysis.csv", index=False)
    print("\n== Power and effect size vs the TRUE fraction of co-timed orthologs")
    print("   (jitter SD 2 h; 60% of genes peaking in one 6 h window)")
    print(t.to_string(index=False))


if __name__ == "__main__":
    main()
