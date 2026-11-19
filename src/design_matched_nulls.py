"""Null calibration under the ACTUAL sampling designs, not one synthetic design.

Two separate validation questions, kept apart:

  (A) our designs  - does the primary null behave correctly under the sampling
      schemes of the cyanobacterial datasets actually analysed here? This is the
      one that licenses the manuscript's statistics.

  (B) published designs - how does a time-point-permutation null behave under
      the sampling schemes reported by Ferrari et al.? This is DESIGN-matched
      only. Their empirical phase distributions and expression characteristics
      are not used, because that requires their data; concentration is swept
      instead. It therefore says what the procedure does under conditions
      resembling their experiments, and says nothing about their actual result.
      A data-matched evaluation belongs to the separate Archaeplastida study.

Inner permutations are set to 1000 for (B) to match the published procedure, and
reported separately from outer trials throughout.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from simcheck import Checkpoint, cell_rng

sys.path.insert(0, str(Path(__file__).parent))
import load
from analyze import OUT
from paths import GEO
from three_nulls_fast import stat_many, draw_marginal, phases_np, mc_ci, PERIOD

# (B) sampling designs as reported by Ferrari et al. 2019
FERRARI = {
    "F: 12 pts, 2 h, 1 day":      np.arange(1, 24, 2.0),
    "F: 7 pts, 4 h (Arab/rice)":  np.arange(0, 25, 4.0),
    "F: 24 pts, 2 h, 2 days":     np.arange(0, 48, 2.0),
}


def our_designs():
    """time vectors actually used by each analysed dataset."""
    out = {}
    for did in ("SE7942_ito_WT", "SYN6803", "TERY", "MED4", "NPUN"):
        _, meta = load.LOADERS[did](GEO)
        out[f"ours: {did}"] = meta["time"].to_numpy(float)
    return out


def one_trial(n, frac, times, rng, n_inner, which, amp=1.0, noise=0.5):
    ca, cb = rng.random() * PERIOD, rng.random() * PERIOD
    pa = draw_marginal(n, frac, 6.0, ca, rng)
    pb = draw_marginal(n, frac, 6.0, cb, rng)       # zero conservation by construction
    t = np.asarray(times, float)
    Ya = amp * np.cos(2 * np.pi * (t[None, :] - pa[:, None]) / PERIOD) + \
        rng.normal(0, noise, (n, len(t)))
    Yb = amp * np.cos(2 * np.pi * (t[None, :] - pb[:, None]) / PERIOD) + \
        rng.normal(0, noise, (n, len(t)))
    fa, fb = phases_np(Ya, t), phases_np(Yb, t)
    obs = stat_many(fa, fb[None, :])[0]
    if which == "marginal":
        null = stat_many(fa, np.array([rng.permutation(fb) for _ in range(n_inner)]))
    else:                                            # time-point permutation
        null = np.empty(n_inner)
        for i in range(n_inner):
            null[i] = stat_many(phases_np(Ya, rng.permutation(t)),
                                phases_np(Yb, rng.permutation(t))[None, :])[0]
    return float((null <= obs).mean())


def run(ck, designs, which, inner, outer, fracs, n, seed):
    for name, times in designs.items():
        for frac in fracs:
            if ck.done((name, which, frac)):
                continue
            rng = cell_rng(seed, (name, which, frac))
            try:
                ps = [one_trial(n, frac, times, rng, inner, which) for _ in range(outer)]
            except Exception as e:
                ck.failed((name, which, frac), e)
                continue
            k = int(sum(p < .05 for p in ps))
            lo, hi = mc_ci(k, outer)
            ck.add(dict(design=name, n_timepoints=len(np.unique(times)),
                        n_samples=len(times), null=which,
                        outer_trials=outer, inner_permutations=inner,
                        frac_in_6h_window=frac, n_orthologs=n,
                        FPR=round(k / outer, 3), FPR_ci95=f"{lo:.3f}-{hi:.3f}"))


def main():
    ck = Checkpoint(OUT, "n16_design_matched_nulls.csv", ["design", "null", "frac_in_6h_window"])
    # (A) our designs, primary null - the calibration the manuscript relies on
    run(ck, our_designs(), "marginal", inner=1000, outer=300,
        fracs=(0.0, 0.6), n=300, seed=31)
    # (B) published designs, time-point permutation at the published 1000 permutations
    run(ck, FERRARI, "timepoint", inner=1000, outer=300,
        fracs=(0.0, 0.4, 0.6, 0.9), n=300, seed=37)
    exp = ([(d, "marginal", f) for d in our_designs() for f in (0.0, 0.6)]
           + [(d, "timepoint", f) for d in FERRARI for f in (0.0, 0.4, 0.6, 0.9)])
    t = ck.finish(sort_by=["null", "design", "frac_in_6h_window"], expected=exp,
                  config=dict(outer_trials=300, inner_permutations=1000,
                              seed_our=31, seed_published=37, n_orthologs=300,
                              seeding="per-cell from (seed, design, null, frac)",
                              designs_ours={k: list(map(float, v))
                                            for k, v in our_designs().items()},
                              designs_published={k: list(map(float, v))
                                                 for k, v in FERRARI.items()},
                              alpha=0.05))
    print("\n== false-positive rate by sampling design (zero conservation by construction)")
    print(t.to_string(index=False))


if __name__ == "__main__":
    main()
