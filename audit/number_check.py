#!/usr/bin/env python3
"""Pre-release number check (RULES.md §5: "the results file wins").

RULES.md names `audit/number_check*.py` as the automated check that re-computes every number in
the text before a release. This project's checking machinery lives in `src/` and runs inside the
build; this script is the single release-gate entry point the rules describe. It does not
re-implement the checks -- duplicating them would create a second source of truth, which is the
thing they exist to prevent. It runs them and reports one verdict.

What it runs:

  src/check_manuscript.py   every statistic in the Abstract, Results and Discussion is re-derived
                            from results/ and compared with what the manuscript prints; any
                            hand-typed statistic fails
  src/validate_semantics.py the claim-level assertions: identity, provenance, claim strength, the
                            freeze conditions, reproducibility
  src/number_provenance.py  every registry key traced to the result file that produced it

What a pass means, and does not mean: the numbers in the text match the committed results, and the
claim-level assertions hold. It does NOT mean the inference is right. A reproducible computation can
reproducibly answer the wrong question -- see CLOSURE-REPORT.md, and §3.9 of the manuscript for the
worked example.

Writes audit/number_check.tsv: every statistic in the manuscript with the file it came from.
Exit status 0 only if every check passed.

usage (from the project directory): python3 audit/number_check.py
"""
import csv
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

CHECKS = [
    ("numbers match results/, none hand-typed", "src/check_manuscript.py"),
    ("claim-level semantic assertions", "src/validate_semantics.py"),
    ("every statistic traced to its source file", "src/number_provenance.py"),
]


def run(label, script):
    print(f"\n=== {label}\n    {script}")
    r = subprocess.run([sys.executable, script], cwd=ROOT)
    print(f"    -> {'PASS' if r.returncode == 0 else 'FAIL'}")
    return r.returncode == 0


def write_table():
    """One row per manuscript statistic: the value, and the file that produced it."""
    src = ROOT / "results" / "NUMBER_PROVENANCE.csv"
    if not src.exists():
        print(f"\n{src.name} not found; run `make validate` first")
        return 0
    with src.open() as fh:
        rows = list(csv.DictReader(fh))
    out = HERE / "number_check.tsv"
    cols = ["key", "value", "source_files", "producer_target", "resolution"]
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for r in sorted(rows, key=lambda r: r["key"]):
            w.writerow([r.get(c, "") for c in cols])
    print(f"\nwrote {out.relative_to(ROOT)}  ({len(rows)} statistics)")
    return len(rows)


def main():
    results = [(label, run(label, script)) for label, script in CHECKS]
    n = write_table()
    failed = [label for label, ok in results if not ok]
    print("\n" + "=" * 62)
    for label, ok in results:
        print(f"  {'PASS' if ok else 'FAIL'}  {label}")
    if failed:
        print(f"\nNUMBER CHECK FAILED: {len(failed)} of {len(results)} checks")
        return 1
    print(f"\nNUMBER CHECK OK: {n} statistics, all re-derived from results/")
    print("A pass means the numbers match the committed results. It does not mean")
    print("the inference is right; CLOSURE-REPORT.md separates those verdicts.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
