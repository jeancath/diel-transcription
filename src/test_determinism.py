"""Interrupt a simulation, resume it, and prove the output is unchanged.

The per-cell seeding in simcheck.py is supposed to make a resumed run
identical to an uninterrupted one. That is a claim about behaviour, so it is
tested by behaviour: run a grid to completion, then run the same grid in a
fresh process, kill it mid-way, restart it, and compare the two final CSVs
cell by cell.

A small grid is used so the test is affordable. The property under test -
whether a cell's value depends on how many cells ran before it - does not
depend on grid size.

Run: python3 src/test_determinism.py
"""
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
WORK = ROOT / "results" / "_determinism"

# A grid small enough to run twice, large enough to interrupt in the middle.
# Negative control: one RNG stream advanced across cells, i.e. the design this
# replaced. A resumed run must then DIFFER, because skipped cells no longer
# consume draws. If the control does not differ, the test is not sensitive and
# its pass means nothing.
RUNNER_SHARED_RNG = """
import sys
sys.path.insert(0, {src!r})
import numpy as np
from simcheck import Checkpoint
from three_nulls_fast import one_trial

ck = Checkpoint({work!r}, "grid.csv", ["n_orthologs", "frac_in_6h_window"])
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
"""

RUNNER = '''
import sys
sys.path.insert(0, {src!r})
import numpy as np
from simcheck import Checkpoint, cell_rng
from three_nulls_fast import one_trial

ck = Checkpoint({work!r}, "grid.csv", ["n_orthologs", "frac_in_6h_window"])
times = np.arange(0, 48, 2.0)
for n in (120, 200):
    for frac in (0.0, 0.3, 0.6, 0.9):
        if ck.done((n, frac)):
            continue
        rng = cell_rng(99, (n, frac))
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
'''


def run(tag, kill_after_cells=None, runner=None):
    """Run the grid. If kill_after_cells is set, SIGKILL once that many cells
    have been written, then return without finishing."""
    script = WORK / f"runner_{tag}.py"
    script.write_text((runner or RUNNER).format(src=str(SRC), work=str(WORK)))
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    p = subprocess.Popen([sys.executable, str(script)], env=env,
                         stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if kill_after_cells is None:
        _, err = p.communicate(timeout=900)
        if p.returncode != 0:
            raise RuntimeError(err.decode()[-2000:])
        return "completed"
    partial = WORK / "grid.partial.csv"
    deadline = time.time() + 900
    while time.time() < deadline:
        if p.poll() is not None:
            raise RuntimeError("process finished before it could be interrupted")
        if partial.exists():
            try:
                if len(pd.read_csv(partial)) >= kill_after_cells:
                    break
            except Exception:
                pass                      # mid-write; try again
        time.sleep(0.2)
    os.kill(p.pid, signal.SIGKILL)
    p.wait()
    return "killed"


def negative_control():
    """The same experiment under the shared-RNG design, which must differ."""
    global WORK
    ctl = WORK / "control"
    keep, WORK = WORK, ctl
    try:
        if ctl.exists():
            shutil.rmtree(ctl)
        ctl.mkdir(parents=True)
        run("ca", runner=RUNNER_SHARED_RNG)
        a = pd.read_csv(ctl / "grid.csv")
        (ctl / "grid.csv").unlink()
        run("cb", kill_after_cells=3, runner=RUNNER_SHARED_RNG)
        run("cb2", runner=RUNNER_SHARED_RNG)
        b = pd.read_csv(ctl / "grid.csv")
    finally:
        WORK = keep
    key = ["n_orthologs", "frac_in_6h_window"]
    a = a.sort_values(key).reset_index(drop=True)
    b = b.sort_values(key).reset_index(drop=True)
    same = all(a[c].equals(b[c]) for c in a.columns)
    n_diff = sum(1 for c in a.columns for i in range(len(a))
                 if a[c].iloc[i] != b[c].iloc[i])
    print(f"  shared-RNG design: resumed run identical? {same} "
          f"({n_diff} differing values)")
    return not same, n_diff


def main():
    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True)

    print("== run A: uninterrupted")
    run("a")
    a = pd.read_csv(WORK / "grid.csv")
    (WORK / "grid.csv").rename(WORK / "grid_uninterrupted.csv")
    print(f"   {len(a)} cells")

    print("== run B: killed mid-grid, then resumed")
    run("b", kill_after_cells=3)
    partial = WORK / "grid.partial.csv"
    n_before = len(pd.read_csv(partial))
    print(f"   killed with {n_before} cells checkpointed")
    # the kill must not have left a corrupt checkpoint
    try:
        pd.read_csv(partial)
        readable = True
    except Exception:
        readable = False
    run("b2")
    b = pd.read_csv(WORK / "grid.csv")
    print(f"   resumed to {len(b)} cells")

    key = ["n_orthologs", "frac_in_6h_window"]
    a = a.sort_values(key).reset_index(drop=True)
    b = b.sort_values(key).reset_index(drop=True)

    fails = []
    if not readable:
        fails.append("checkpoint left unreadable by the kill")
    if n_before < 1 or n_before >= len(a):
        fails.append(f"interruption was not mid-grid ({n_before} of {len(a)})")
    if list(a.columns) != list(b.columns) or len(a) != len(b):
        fails.append("shape differs between runs")
    else:
        diff = []
        for c in a.columns:
            if not a[c].equals(b[c]):
                for i in range(len(a)):
                    if a[c].iloc[i] != b[c].iloc[i]:
                        diff.append(f"{c} row {i}: {a[c].iloc[i]} vs {b[c].iloc[i]}")
        if diff:
            fails.append("cell values differ: " + "; ".join(diff[:6]))

    print()
    print(f"  cells compared:            {len(a)}")
    print(f"  cells written before kill: {n_before}")
    print(f"  cells recomputed after:    {len(a) - n_before}")
    print(f"  checkpoint readable after SIGKILL: {readable}")
    print(f"  identical, every cell and column:  {not fails}")
    a.to_csv(WORK / "grid_uninterrupted.csv", index=False)
    b.to_csv(WORK / "grid_resumed.csv", index=False)
    print(f"\n  evidence: {WORK.relative_to(ROOT)}/grid_uninterrupted.csv, grid_resumed.csv")

    print("\n== negative control: the old shared-RNG design")
    differs, n_diff = negative_control()
    if not differs:
        fails.append("negative control did not differ - the test is insensitive")
    else:
        print(f"  control differs in {n_diff} values, so the comparison "
              f"above is capable of failing")

    import json
    rep = {
        "cells_compared": int(len(a)),
        "cells_written_before_kill": int(n_before),
        "cells_recomputed_after_resume": int(len(a) - n_before),
        "checkpoint_readable_after_sigkill": bool(readable),
        "resumed_equals_uninterrupted": not fails,
        "negative_control_shared_rng_differs": bool(differs),
        "negative_control_differing_values": int(n_diff),
        "grid": "three_nulls_fast.one_trial, n in (120, 200) x "
                "concentration in (0.0, 0.3, 0.6, 0.9), 25 outer x 40 inner",
        "seed": 99,
        "failures": fails,
    }
    (ROOT / "results" / "determinism_report.json").write_text(
        json.dumps(rep, indent=2) + "\n")
    print(f"  report: results/determinism_report.json")

    if fails:
        print("\nDETERMINISM TEST FAILED")
        for f in fails:
            print(f"  - {f}")
        return 1
    print("\nDETERMINISM OK: resumed run is bit-identical to the uninterrupted run")
    return 0


if __name__ == "__main__":
    sys.exit(main())
