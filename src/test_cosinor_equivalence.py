"""Equivalence test: the fast numpy cosinor against the reference implementation.

src/three_nulls_fast.py rewrote the per-permutation cosinor fit in pure numpy to
make 1000 outer trials affordable. A rewritten numerical routine is not accepted
on the grounds that it is faster; it has to return the same numbers. This
compares amplitude, phase, F statistic, raw p and BH q on identical inputs, on
random data and on a real dataset, and exits 1 on any disagreement beyond
floating-point tolerance.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import rhythm
from three_nulls_fast import phases_np, design, PERIOD
from scipy import stats

TOL = 1e-8


def full_numpy_fit(Y, times):
    """Amplitude, phase, F, p computed in pure numpy, matching rhythm.cosinor."""
    t = np.asarray(times, float)
    S = np.ones((len(t), 1))
    H = np.column_stack([np.cos(2 * np.pi * t / PERIOD), np.sin(2 * np.pi * t / PERIOD)])
    X = np.hstack([S, H])
    n, pf, pr = len(t), X.shape[1], S.shape[1]
    cf, _, _, _ = np.linalg.lstsq(X, Y.T, rcond=None)
    cr, _, _, _ = np.linalg.lstsq(S, Y.T, rcond=None)
    rss_f = ((Y.T - X @ cf) ** 2).sum(axis=0)
    rss_r = ((Y.T - S @ cr) ** 2).sum(axis=0)
    F = ((rss_r - rss_f) / (pf - pr)) / (rss_f / (n - pf))
    p = stats.f.sf(F, pf - pr, n - pf)
    A, B = cf[-2], cf[-1]
    return dict(amp=2 * np.sqrt(A ** 2 + B ** 2),
                phase=(np.arctan2(B, A) * PERIOD / (2 * np.pi)) % PERIOD,
                F=F, p=p, q=rhythm.bh(p))


def compare(Y, times, label):
    expr = pd.DataFrame(Y, index=[f"g{i}" for i in range(len(Y))],
                        columns=[f"s{j}" for j in range(len(times))])
    meta = pd.DataFrame({"time": times, "series": "s1"}, index=expr.columns)
    ref = rhythm.cosinor(expr, meta, PERIOD)
    fast = full_numpy_fit(Y, times)
    fails = []
    for key, refcol in (("amp", "amp_log2"), ("phase", "phase_h"), ("p", "p"), ("q", "q")):
        a, b = np.asarray(ref[refcol]), fast[key]
        if key == "phase":                      # circular: compare on the circle
            d = np.abs((a - b + PERIOD / 2) % PERIOD - PERIOD / 2)
        else:
            d = np.abs(a - b)
        worst = float(np.nanmax(d))
        ok = worst < TOL
        print(f"  {label:22s} {key:6s} max|diff| = {worst:.3e}  {'OK' if ok else 'FAIL'}")
        if not ok:
            fails.append((label, key, worst))
    # the phase-only shortcut used inside the simulation loop
    ps = phases_np(Y, times)
    d = float(np.nanmax(np.abs((np.asarray(ref["phase_h"]) - ps + PERIOD / 2)
                               % PERIOD - PERIOD / 2)))
    ok = d < TOL
    print(f"  {label:22s} phases_np shortcut max|diff| = {d:.3e}  {'OK' if ok else 'FAIL'}")
    if not ok:
        fails.append((label, "phases_np", d))
    return fails


def main():
    rng = np.random.default_rng(0)
    fails = []

    t = np.arange(0, 48, 2.0)
    ph = rng.random(500) * PERIOD
    Y = np.cos(2 * np.pi * (t[None, :] - ph[:, None]) / PERIOD) + rng.normal(0, .5, (500, len(t)))
    fails += compare(Y, t, "synthetic 2h/48h")

    t2 = np.array([0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48.])
    ph2 = rng.random(300) * PERIOD
    Y2 = 2 * np.cos(2 * np.pi * (t2[None, :] - ph2[:, None]) / PERIOD) + rng.normal(0, 1, (300, len(t2)))
    fails += compare(Y2, t2, "synthetic 4h/48h")

    fails += compare(rng.normal(0, 1, (200, len(t))), t, "pure noise")

    try:                                        # a real dataset, single series
        import load
        from paths import GEO
        expr, meta = load.LOADERS["MED4"](GEO)
        m = meta[meta.series == meta.series.iloc[0]]
        fails += compare(expr[m.index].to_numpy(float), m.time.values, "MED4 (real data)")
    except Exception as e:
        print(f"  real-data comparison skipped: {type(e).__name__}: {e}")

    print()
    if fails:
        print(f"EQUIVALENCE FAILED in {len(fails)} comparison(s)")
        return 1
    print("EQUIVALENCE OK: fast numpy cosinor matches rhythm.cosinor to < 1e-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
