"""A provenance record for every number the manuscript prints.

For each registry key: where it appears in the manuscript, which result file it
came from, which script and make target produced that file, and the test that
proves the printed number actually depends on it.

The mapping is derived by static analysis of manuscript_numbers.collect(), not
maintained by hand. A hand-kept ledger is the thing that goes stale: the §3.4
power prose quoted a 40-trial pilot for a whole version while the table beside
it reported the 400-trial surface, and nothing noticed.

Method. Walk the registry source in order. Any local name assigned from
_csv("f.csv") or pd.read_csv(...) carries those filenames. For every
`v["key"] = expr`, the names used in expr are resolved transitively back to the
files they came from. Keys built from f-strings (v[f"within_{rule}"]) are
matched to real keys by pattern. Keys that resolve to no file are reported as
undeclared and fail the run unless they are on CONSTANTS, which holds values
this project chooses rather than measures.

The dependency test is a perturbation, the same device the renderer audit uses:
change the key's value, re-render, and require the output to change. That is
what distinguishes a documented number from a number that merely sits in a
dictionary.
"""
import ast
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import manuscript_numbers as MN
import render_manuscript as R
from analyze import OUT

SRC = Path(__file__).resolve().parent
ROOT = SRC.parent

# Registry values this project sets rather than measures. Each needs a reason.
CONSTANTS = {
    "alpha": "nominal significance level, chosen",
    "qcut": "rhythm-call FDR threshold, chosen",
    "period": "assumed circadian period",
}


def _files_in(node):
    """Result-file names referenced anywhere inside a node.

    Both spellings are picked up: _csv("f.csv") / read_csv(...) calls, and a
    bare filename constant such as OUT / "n18_x.csv", which is often assigned
    to a local on one line and read on another.
    """
    out = []
    for n in ast.walk(node):
        if isinstance(n, ast.Constant) and isinstance(n.value, str) and \
                n.value.endswith((".csv", ".csv.gz", ".json")) and \
                not n.value.startswith((".", "_")):
            # a leading "." or "_" means an f-string fragment such as
            # f"n18_stratified_closure_{lvl}_last.csv", not a whole filename
            out.append(n.value)
    return out


def _names_in(node):
    return {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}


def _key_pattern(node):
    """The literal key, or a regex for an f-string key."""
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value, None
    if isinstance(node, ast.JoinedStr):
        pat = ""
        for part in node.values:
            if isinstance(part, ast.Constant):
                pat += re.escape(str(part.value))
            else:
                pat += r"[A-Za-z0-9_.]+"
        return None, re.compile("^" + pat + "$")
    return None, None


def trace_registry(func_name="collect"):
    """key (or key pattern) -> set of result files it was computed from."""
    tree = ast.parse((SRC / "manuscript_numbers.py").read_text())
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == func_name)
    owner = {}       # local name -> files
    literal, dynamic = {}, []
    # Source order, not ast.walk order. ast.walk is breadth-first, so a
    # statement inside an `if` block is visited after the top-level statements
    # that consume its result, and the name it binds looks untraceable.
    # Keys written as v[tag] inside `for col, tag in (...)` are resolved to the
    # string constants the loop iterates over; keys written as v[f"..."] to a
    # pattern. Both appear in the registry and neither is a literal.
    loopvars = {}
    for n in ast.walk(fn):
        if isinstance(n, ast.For):
            names = {x.id for x in ast.walk(n.target) if isinstance(x, ast.Name)}
            lits = [c.value for c in ast.walk(n.iter)
                    if isinstance(c, ast.Constant) and isinstance(c.value, str)]
            for nm in names:
                loopvars.setdefault(nm, set()).update(lits)

    # A list or set built up with .append()/.update() carries the provenance of
    # everything appended to it; those are Expr statements, not assignments.
    mutations_ = sorted((n for n in ast.walk(fn) if isinstance(n, ast.Expr)
                         and isinstance(n.value, ast.Call)
                         and isinstance(n.value.func, ast.Attribute)
                         and n.value.func.attr in ("append", "extend", "update", "add")),
                        key=lambda n: (n.lineno, n.col_offset))
    assigns = sorted((n for n in ast.walk(fn)
                      if isinstance(n, (ast.Assign, ast.AugAssign))) ,
                     key=lambda n: (n.lineno, n.col_offset))
    # `for _rule, _g in ep.groupby("rule")` binds _g to a slice of ep, so the
    # loop target inherits the iterable's provenance.
    fors = sorted((n for n in ast.walk(fn) if isinstance(n, ast.For)),
                  key=lambda n: (n.lineno, n.col_offset))
    stmts = sorted(assigns + mutations_ + fors,
                   key=lambda n: (n.lineno, n.col_offset))
    for node in stmts:
        if isinstance(node, ast.For):
            files = set(_files_in(node.iter))
            for nm in _names_in(node.iter):
                files |= owner.get(nm, set())
            for x in ast.walk(node.target):
                if isinstance(x, ast.Name):
                    owner.setdefault(x.id, set()).update(files)
            continue
        if isinstance(node, ast.Expr):
            holder = getattr(node.value.func.value, "id", None)
            if holder:
                files = set(_files_in(node.value))
                for nm in _names_in(node.value):
                    files |= owner.get(nm, set())
                owner.setdefault(holder, set()).update(files)
            continue
        targets = node.targets if isinstance(node, ast.Assign) else [node.target]
        files = set(_files_in(node.value))
        for nm in _names_in(node.value):
            files |= owner.get(nm, set())
        for t in targets:
            # v["key"] = ... / v[f"..."] = ...
            if isinstance(t, ast.Subscript) and getattr(t.value, "id", "") == "v":
                k, pat = _key_pattern(t.slice)
                if k is not None:
                    literal.setdefault(k, set()).update(files)
                elif pat is not None:
                    dynamic.append((pat, files))
                elif isinstance(t.slice, ast.Name):
                    for lit in loopvars.get(t.slice.id, ()):
                        literal.setdefault(lit, set()).update(files)
            elif isinstance(t, ast.Name):
                owner[t.id] = files
            elif isinstance(t, ast.Tuple):
                for el in t.elts:
                    if isinstance(el, ast.Name):
                        owner[el.id] = files
                    elif isinstance(el, ast.Subscript) and \
                            getattr(el.value, "id", "") == "v":
                        k, pat = _key_pattern(el.slice)
                        if k is not None:
                            literal.setdefault(k, set()).update(files)
                        elif pat is not None:
                            dynamic.append((pat, files))
    return literal, dynamic


def main():
    values = MN.collect()
    tables = MN.tables(values)
    tpl = R.TEMPLATE.read_text()
    man_lines = R.TARGET.read_text().splitlines()
    literal, dynamic = trace_registry()

    # where each placeholder sits in the template, by line
    loc = {}
    for i, line in enumerate(tpl.splitlines(), 1):
        for k in re.findall(R.PLACEHOLDER, line):
            loc.setdefault(k, []).append(i)

    baseline = R.render(tpl, values=values, tables=tables)

    def perturb(x):
        if isinstance(x, bool):
            return not x
        if isinstance(x, (int, float)):
            return type(x)(x + 7919)
        return str(x) + "_MUT7919"

    rows, undeclared = [], []
    for k in sorted(values):
        files = set(literal.get(k, set()))
        how = "direct" if files else ""
        if not files:
            for pat, f in dynamic:
                if pat.match(k):
                    files |= f
                    how = "pattern"
        if not files and k not in CONSTANTS:
            undeclared.append(k)
            how = "UNDECLARED"
        elif not files:
            how = "constant"

        mut = dict(values)
        mut[k] = perturb(values[k])
        try:
            changed = R.render(tpl, values=mut, tables=MN.tables(mut)) != baseline
        except Exception:
            changed = True       # a formatting error means the key is read

        prod = sorted({MN.producer(f) or "?" for f in files}) if files else []
        rows.append(dict(
            key=k, value=str(values[k]),
            template_lines=";".join(map(str, loc.get(k, []))),
            n_occurrences=len(loc.get(k, [])),
            appears_as=("placeholder" if loc.get(k)
                        else "table generator input"),
            source_files=";".join(sorted(files)),
            resolution=how,
            producer_target=";".join(p for p in prod if p),
            output_depends_on_key=changed,
            note=CONSTANTS.get(k, "")))

    t = pd.DataFrame(rows)
    t.to_csv(OUT / "NUMBER_PROVENANCE.csv", index=False)
    print(f"registry keys:                 {len(t)}")
    print(f"traced to a result file:       {int((t.source_files != '').sum())}")
    print(f"  directly                     {int((t.resolution == 'direct').sum())}")
    print(f"  by f-string key pattern      {int((t.resolution == 'pattern').sum())}")
    print(f"declared constants:            {int((t.resolution == 'constant').sum())}")
    print(f"printed number depends on key: {int(t.output_depends_on_key.sum())}")
    print(f"written {OUT / 'NUMBER_PROVENANCE.csv'}")
    bad = t[~t.output_depends_on_key]
    if len(bad):
        print("\nFAIL  keys the rendered manuscript does not depend on: "
              f"{bad.key.tolist()}")
    if undeclared:
        print(f"\nFAIL  {len(undeclared)} key(s) trace to no result file and are "
              f"not declared constants:\n  {undeclared}")
    return 1 if (undeclared or len(bad)) else 0


if __name__ == "__main__":
    sys.exit(main())
