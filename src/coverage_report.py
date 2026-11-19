"""Map semantic assertions to the mutation tests that prove they fire.

29 assertions and 10 mutation tests are two different quantities, and quoting
them together invites the reading that each assertion has a test. Most do not.
This script derives both lists from the source - the `check(...)` calls in
validate_semantics.py and the `case(...)` registrations in
test_semantic_mutations.py - matches each case to the assertion its expected
substring selects, and writes the mapping to results/ASSERTION_COVERAGE.csv.

A mutation test counts as covering an assertion only if the validator fails on
that named assertion. test_semantic_mutations.py enforces that at run time by
requiring the named assertion to be the one that failed, so a generic crash
cannot be mistaken for coverage; this script reports which assertions have such
a test and which rest on review alone.
"""
import ast
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from analyze import OUT

SRC = Path(__file__).resolve().parent


# Every assertion site is classified, so that "14 sites rest on review" can be
# read as a risk statement rather than a count. The categories and the action
# each implies were set by review; a site matching no rule is reported as
# unclassified rather than silently defaulted, because a new critical check
# must not land in a lenient bucket by accident.
CATEGORIES = (
    ("kaib and kaic are different genes",      "critical: biological identity"),
    ("clock-gene label matches",               "critical: biological identity"),
    ("every analysed arm appears",             "critical: biological identity"),
    ("gene symbol named in 3.9",               "critical: biological identity"),
    ("accession",                              "critical: biological identity"),
    ("each downloaded matrix declares",        "critical: biological identity"),
    ("fitted phases lie in",                   "critical: statistical claim"),
    ("q values lie in",                        "critical: statistical claim"),
    ("kaiabc-null control",                    "critical: statistical claim"),
    ("reported rbh range",                     "critical: statistical claim"),
    ("claims stronger than",                   "critical: statistical claim"),
    ("distinguishes diel expression",          "critical: statistical claim"),
    ("design-matched only",                    "critical: statistical claim"),
    ("no claim wider than",                    "critical: statistical claim"),
    ("qualifications its claims require",      "critical: statistical claim"),
    ("brite assignment rule",                  "critical: statistical claim"),
    ("neither brite rule was chosen",          "critical: statistical claim"),
    ("enrichment is not established",          "critical: statistical claim"),
    ("sensitivity, not absence",               "critical: statistical claim"),
    ("endogenous circadian scope",             "critical: statistical claim"),
    ("abstract states the study is exploratory", "critical: statistical claim"),
    ("abstract states the diel versus",         "critical: statistical claim"),
    ("no p value or fold enrichment",          "critical: statistical claim"),
    ("series count is derived",                "structural / documentation"),
    ("arms outnumber series",                  "structural / documentation"),
    ("unresolved template placeholder",         "structural / documentation"),
    ("declared by fetch_data",                 "structural / documentation"),
    ("symlink outside the repository",         "structural / documentation"),
    ("producing script",                       "structural / documentation"),
    ("producing target is reachable",          "structural / documentation"),
    ("definition follows the __main__",        "structural / documentation"),
    ("symlink portability check could run",    "low impact: meta-check"),
)
ACTION = {
    "critical: biological identity": "add a targeted mutation test",
    "critical: statistical claim": "add a targeted mutation test",
    "structural / documentation": "documented review may suffice",
    "low impact: meta-check": "review accepted",
    "UNCLASSIFIED": "classify, then decide",
}


def categorise(name):
    low = name.lower()
    for key, cat in CATEGORIES:
        if key in low:
            return cat
    return "UNCLASSIFIED"


def _const_str(node):
    return node.value if isinstance(node, ast.Constant) and \
        isinstance(node.value, str) else None


def assertions():
    """Every assertion name, in the order the validator checks them."""
    tree = ast.parse((SRC / "validate_semantics.py").read_text())
    out = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) \
                and n.func.id == "check" and n.args:
            s = _const_str(n.args[0])
            if s:
                out.append((n.lineno, s))
    return [s for _, s in sorted(out)]


def mutations():
    """Each mutation case: label, mutator function, expected assertion substring."""
    tree = ast.parse((SRC / "test_semantic_mutations.py").read_text())
    out = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) \
                and n.func.id == "case" and len(n.args) == 3:
            label, fn, expect = n.args
            out.append(dict(mutation=_const_str(label),
                            mutator=getattr(fn, "id", "?"),
                            expects=_const_str(expect)))
    return [m for m in out if m["mutation"] and m["expects"]]


def main():
    asserts, muts = assertions(), mutations()
    rows = []
    for a in asserts:
        hit = [m for m in muts if m["expects"] and m["expects"].lower() in a.lower()]
        cat = categorise(a)
        rows.append(dict(assertion=a, category=cat,
                         validation=("targeted mutation test" if hit
                                     else "review only"),
                         covered_by_mutation="; ".join(m["mutation"] for m in hit),
                         injected_error="; ".join(m["mutator"] for m in hit),
                         covered=bool(hit),
                         action="" if hit else ACTION[cat]))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "ASSERTION_COVERAGE.csv", index=False)

    orphan = [m["mutation"] for m in muts
              if not any(m["expects"].lower() in a.lower() for a in asserts)]
    print(f"semantic assertions:         {len(asserts)}")
    print(f"mutation tests:              {len(muts)}")
    print(f"assertions with a mutation test: {int(t.covered.sum())}")
    print(f"assertions resting on review:    {int((~t.covered).sum())}")
    print("\n== by category: proven by a targeted mutation vs resting on review")
    print(f"{'category':32s} {'mutation':>9s} {'review':>7s}")
    for cat, g in t.groupby("category"):
        print(f"{cat:32s} {int(g.covered.sum()):9d} {int((~g.covered).sum()):7d}")
    crit = t[(~t.covered) & t.category.str.startswith("critical")]
    if len(crit):
        print(f"\n!! {len(crit)} CRITICAL site(s) rest on review alone and are NOT "
              f"experimentally validated:")
        for r in crit.itertuples():
            print(f"   [{r.category}] {r.assertion}")
        print("   action: " + ACTION["critical: statistical claim"]
              + " (recorded as follow-up work, not done for v0.3)")
    un = t[t.category == "UNCLASSIFIED"]
    if len(un):
        print(f"\nFAIL  unclassified assertion site(s): {un.assertion.tolist()}")
        return 1
    print("\n== assertions proven to fire by a targeted mutation")
    for r in t[t.covered].itertuples():
        print(f"  {r.assertion}\n      <- {r.covered_by_mutation}  [{r.injected_error}]")
    if orphan:
        print("\n!! mutation tests matching no assertion name: " + "; ".join(orphan))
        return 1
    print(f"\nwritten {OUT / 'ASSERTION_COVERAGE.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
