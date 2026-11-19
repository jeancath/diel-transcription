"""Compare two builds file by file and classify every difference.

Usage: python3 src/compare_builds.py OLD_RESULTS_DIR NEW_RESULTS_DIR

A build that "succeeds" is not evidence of reproducibility; identical outputs
are. Each file is compared by content and every difference is classified:

  identical            byte-for-byte, or numerically equal within tolerance
  expected: <reason>   differs for a declared reason (a timestamp, the commit,
                       wall-clock timing, or a live external source)
  UNEXPECTED           anything else - these are the ones that matter

Numeric CSV columns are compared with a tolerance so that a difference in the
last printed digit is not reported as a change in the science, while any
larger difference is.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Fields and files that are expected to differ between two builds, with the
# reason. Anything not listed here must be identical.
VOLATILE_JSON_KEYS = {
    "written_utc", "wall_seconds", "git_commit", "commit", "describe",
    "tree_clean", "dirty_paths", "data_dir", "platform", "machine",
}
LIVE_SOURCES = {
    "lit_query_log.csv": "PubMed record counts move as records are indexed",
    "lit_screen.csv": "PubMed result sets move as records are indexed",
    "lit_ferrari_citing_screen.csv": "the citation graph grows over time",
}
TOL = 1e-9

# Files whose BYTES differ between identical builds for a reason outside the
# analysis. gzip writes the current time into its header, so a .csv.gz is not
# byte-reproducible even when every value in it is identical; the content is
# compared by reading through the compression, and only the hash recorded in
# PROVENANCE.json sees the difference. Fixing it means passing mtime=0 to the
# writer; it is declared rather than fixed, because the content is verified
# and re-running the build to change a header is not worth it.
NON_BYTE_REPRODUCIBLE = {
    "results/rhythms_per_gene.csv.gz": "gzip stores the current time in its header",
}


# Differences against the PRE-REBUILD BACKUP that are accounted for, each with
# the change that caused it. Only applied with --baseline, and only to the file
# named: anything not listed here is still UNEXPECTED. The point is that a
# difference has to be explained by name, not waved through as "probably fine".
BASELINE_DELTAS = {
    # artifacts of validation and provenance work added after the backup
    "PROVENANCE.json": "new: src/provenance.py added (input/output checksums)",
    "NUMBER_PROVENANCE.csv": "new: src/number_provenance.py added",
    "ASSERTION_COVERAGE.csv": "new: src/coverage_report.py added",
    "determinism_report.json": "new: src/test_determinism.py added (make tests-slow)",
    "n3_rbh_params.json": "new: rbh_check.write_params_record added",
    "n13_three_nulls_1000.manifest.json": "new: simcheck checkpointing added",
    "n14_power_surface.manifest.json": "new: simcheck checkpointing added",
    "n16_design_matched_nulls.manifest.json": "new: simcheck checkpointing added",
    "n17_power_small_fractions.manifest.json": "new: simcheck checkpointing added",
    # analyses added after the backup
    "n5_core_sd_null.csv": "new: sd_null wired into conserved_core (3.9 control)",
    "n18_stratified_closure_A_last.csv": "new: multi='last' BRITE sensitivity run",
    "n18_stratified_closure_B_last.csv": "new: multi='last' BRITE sensitivity run",
    "n19_brite_multiassign.csv": "new: multiassign_report, previously unreachable",
    # schema changes
    "n18_stratified_closure_A.csv": "schema: multi_rule and n_permutations recorded",
    "n18_stratified_closure_B.csv": "schema: multi_rule and n_permutations recorded",
    "n1_two_nulls.csv": "schema: two columns dropped that no committed code wrote "
                        "(effect_vs_marginal_pct, effect_vs_time_null_pct came from "
                        "an uncommitted interactive variant); nothing reads them",
    # deliberate removal
    "n9_three_nulls.csv": "removed: the 50-trial pilot and three_nulls.py deleted; "
                          "make clean deleted it and make all never restored it",
    # recomputation under a changed RNG scheme
    "n17_power_small_fractions.csv":
        "reseeded: per-cell seeding replaced one advancing RNG stream, so the "
        "cells are independent Monte Carlo draws (18 cells, mean diff +0.006, "
        "max 0.070, against SE ~0.032 for two estimates at 400 trials)",
    # registry keys the manuscript prints
    "numbers.json": "registry: 3.9 inferential keys removed, candidate-set and "
                    "7-point-design keys added",
}


def read_csv(p):
    try:
        return pd.read_csv(p)
    except Exception:
        return None


def compare_csv(a, b):
    """Return (same, note). Numeric columns within TOL count as same."""
    ta, tb = read_csv(a), read_csv(b)
    if ta is None or tb is None:
        return a.read_bytes() == b.read_bytes(), "unparseable; compared as bytes"
    if list(ta.columns) != list(tb.columns):
        return False, f"columns differ: {set(ta.columns) ^ set(tb.columns)}"
    if len(ta) != len(tb):
        return False, f"row count {len(ta)} vs {len(tb)}"
    diffs = []
    for c in ta.columns:
        xa, xb = ta[c], tb[c]
        if pd.api.types.is_numeric_dtype(xa) and pd.api.types.is_numeric_dtype(xb):
            d = np.abs(xa.to_numpy(dtype=float) - xb.to_numpy(dtype=float))
            d = np.where(np.isnan(d), 0.0, d)
            if d.max() > TOL:
                i = int(np.argmax(d))
                diffs.append(f"{c}: max |diff| {d.max():.3g} (row {i}: "
                             f"{xa.iloc[i]} vs {xb.iloc[i]})")
        else:
            # NaN != NaN, so fill with a sentinel before comparing as text;
            # otherwise every missing annotation reads as a difference
            sa_ = xa.where(xa.notna(), "\x00NA").astype(str)
            sb_ = xb.where(xb.notna(), "\x00NA").astype(str)
            ne = (sa_ != sb_)
            if ne.any():
                i = int(np.argmax(ne.to_numpy()))
                diffs.append(f"{c}: {int(ne.sum())} rows differ (row {i}: "
                             f"{xa.iloc[i]!r} vs {xb.iloc[i]!r})")
    return (not diffs), "; ".join(diffs[:3])


def compare_provenance(a, b):
    """PROVENANCE.json is the integrity record, so compare it on its meaning.

    Input checksums must match exactly: a difference there means the two builds
    did not read the same raw data, which no later agreement could excuse.
    Output checksums may differ only for files that are expected to differ -
    the simulation manifests, which record their own wall time, and files that
    are not byte-reproducible for a declared reason. Anything else is a real
    difference in a result, even though the file containing it compared equal.
    """
    import json as _json
    x, y = _json.loads(a.read_text()), _json.loads(b.read_text())
    ia = {k: v.get("sha256") for k, v in x.get("inputs", {}).items()}
    ib = {k: v.get("sha256") for k, v in y.get("inputs", {}).items()}
    if ia != ib:
        bad = sorted(k for k in set(ia) | set(ib) if ia.get(k) != ib.get(k))
        return False, f"INPUT checksums differ: {bad[:6]}"
    oa = {k: v.get("sha256") for k, v in x.get("outputs", {}).items()}
    ob = {k: v.get("sha256") for k, v in y.get("outputs", {}).items()}
    diff = sorted(k for k in set(oa) | set(ob) if oa.get(k) != ob.get(k))
    unexplained = [k for k in diff
                   if not k.endswith(".manifest.json")
                   and k not in NON_BYTE_REPRODUCIBLE]
    if unexplained:
        return False, f"output checksums differ with no declared reason: {unexplained}"
    why = ", ".join(f"{k.split('/')[-1]} ({NON_BYTE_REPRODUCIBLE[k]})"
                    for k in diff if k in NON_BYTE_REPRODUCIBLE)
    n_man = sum(1 for k in diff if k.endswith(".manifest.json"))
    return True, (f"inputs identical; {len(oa)} outputs, "
                  f"{n_man} manifest(s) differ by recorded wall time"
                  + (f"; {why}" if why else ""))


def compare_json(a, b):
    try:
        ja, jb = json.loads(a.read_text()), json.loads(b.read_text())
    except Exception:
        return a.read_bytes() == b.read_bytes(), "unparseable"

    def strip(o):
        if isinstance(o, dict):
            return {k: strip(v) for k, v in o.items() if k not in VOLATILE_JSON_KEYS}
        if isinstance(o, list):
            return [strip(x) for x in o]
        return o
    sa, sb = strip(ja), strip(jb)
    if sa == sb:
        return True, "equal ignoring timestamps, commit and timing"
    keys = sorted(set(sa if isinstance(sa, dict) else {}) |
                  set(sb if isinstance(sb, dict) else {}))
    bad = [k for k in keys if isinstance(sa, dict) and sa.get(k) != sb.get(k)]
    return False, f"differing keys: {bad[:6]}"


def main(old, new, baseline=False):
    old, new = Path(old), Path(new)
    deltas = BASELINE_DELTAS if baseline else {}
    names = sorted({f.name for f in old.iterdir() if f.is_file()} |
                   {f.name for f in new.iterdir() if f.is_file()})
    rows = []
    for n in names:
        a, b = old / n, new / n
        if n.endswith((".partial.csv", ".log")):
            continue          # transient checkpoints and build logs
        if not a.exists():
            rows.append((n, f"accounted: {deltas[n]}" if n in deltas else "UNEXPECTED",
                         "present only in the new build"))
            continue
        if not b.exists():
            rows.append((n, f"accounted: {deltas[n]}" if n in deltas else "UNEXPECTED",
                         "missing from the new build"))
            continue
        if n.endswith((".csv", ".csv.gz")):
            same, note = compare_csv(a, b)
        elif n == "PROVENANCE.json":
            same, note = compare_provenance(a, b)
        elif n.endswith(".json"):
            same, note = compare_json(a, b)
        else:
            same, note = a.read_bytes() == b.read_bytes(), ""
        if same and n == "PROVENANCE.json":
            rows.append((n, "identical", note))
        elif same:
            rows.append((n, "identical", note))
        elif n in LIVE_SOURCES:
            rows.append((n, f"expected: {LIVE_SOURCES[n]}", note))
        elif n.endswith(".manifest.json") or n == "PROVENANCE.json":
            rows.append((n, "expected: records this build's own metadata", note))
        elif n in deltas:
            rows.append((n, f"accounted: {deltas[n]}", note))
        else:
            rows.append((n, "UNEXPECTED", note))

    w = max(len(r[0]) for r in rows) + 2
    print(f"{'file':{w}} verdict")
    print("-" * (w + 60))
    for n, verdict, note in rows:
        print(f"{n:{w}} {verdict}" + (f"\n{'':{w}}   {note}" if note and
                                      verdict == "UNEXPECTED" else ""))
    bad = [r for r in rows if r[1] == "UNEXPECTED"]
    ident = sum(1 for r in rows if r[1] == "identical")
    acc = sum(1 for r in rows if r[1].startswith("accounted:"))
    print(f"\n{ident} identical, {len(rows) - ident - len(bad) - acc} expected to "
          f"differ, {acc} accounted for by a declared change, {len(bad)} unexpected")
    if acc:
        print("\nACCOUNTED-FOR DIFFERENCES (each named, with the change that caused it):")
        for n, verdict, _ in rows:
            if verdict.startswith("accounted:"):
                print(f"  {n}\n      {verdict.split(': ', 1)[1]}")
    unused = sorted(set(deltas) - {r[0] for r in rows})
    if unused:
        print(f"\nFAIL  declared differences that did not occur: {unused}")
        return 1
    if bad:
        print("\nUNEXPECTED DIFFERENCES:")
        for n, _, note in bad:
            print(f"  {n}: {note}")
        return 1
    print("\nBUILDS AGREE: every difference is accounted for")
    return 0


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 2:
        print(__doc__)
        sys.exit(2)
    # --baseline allows ONLY the differences declared in BASELINE_DELTAS, each
    # by name; without it no difference is allowed except the volatile metadata
    # and live sources above, which is the right setting for build-vs-build.
    sys.exit(main(args[0], args[1], baseline="--baseline" in sys.argv))
