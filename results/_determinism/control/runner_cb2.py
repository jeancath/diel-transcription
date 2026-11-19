
import sys
sys.path.insert(0, '<project>/src')
import numpy as np
from simcheck import Checkpoint
from three_nulls_fast import one_trial

ck = Checkpoint('<project>/results/_determinism/control', "grid.csv", ["n_orthologs", "frac_in_6h_window"])
times = np.arange(0, 48, 2.0)
rng = np.random.default_rng(99)          # ONE stream for the whole grid
for n in (120, 200):
    for frac in (0.0, 0.3, 0.6, 0.9):
        if ck.done((n, frac)):
            continue
        res = [one_trial(n, frac, 6.0, times, rng, 40) for _ in range(25)]
        pu, pf, pm, mu, mf, mm = map(np.array, zip(*res))
        ck.add(dict(n_orthologs=n, frac_in_6h_window=frac,
                    fpr_uniform=float((pu < .05).mean()),
                    fpr_timepoint=float((pf < .05).mean()),
                    fpr_marginal=float((pm < .05).mean()),
                    mean_uniform=float(mu.mean()),
                    mean_timepoint=float(mf.mean()),
                    mean_marginal=float(mm.mean())))
ck.finish(sort_by=["n_orthologs", "frac_in_6h_window"])
