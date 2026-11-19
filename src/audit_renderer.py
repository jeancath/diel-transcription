"""Independent audit of the manuscript renderer, with mutation tests.

Why this exists. render_manuscript.py used one regular expression both to find
placeholders and to report the ones it could not resolve. The class was
lowercase-only, so a key containing a capital letter was invisible to both
roles: not substituted, not reported. Twenty-one lines of the manuscript
shipped as literal text and every check passed.

The lesson is that a check sharing an implementation with the thing it checks
proves nothing. This file deliberately re-derives everything from the files
themselves - it does not import the renderer's pattern - and then mutates the
template to confirm each check actually fires.

Run: python3 src/audit_renderer.py
"""
import json
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import manuscript_numbers as MN
import render_manuscript as R

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "MANUSCRIPT.template.md"
TARGET = ROOT / "MANUSCRIPT.md"

# Deliberately permissive: anything that looks remotely like a placeholder.
# Written independently of render_manuscript.PLACEHOLDER so that a defect in
# that pattern cannot hide from this one.
ANY_BRACED = re.compile(r"\{\{(.*?)\}\}", re.S)
WELL_FORMED = re.compile(r"^[A-Za-z0-9_]+(?::[A-Za-z0-9_]+)?$")

FAILS = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  — {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)
    return ok


def audit():
    tpl = TEMPLATE.read_text()
    values = MN.collect()
    tables = MN.tables()
    tokens = [m.group(1) for m in ANY_BRACED.finditer(tpl)]

    # 1. every braced token is a well-formed placeholder
    bad = sorted({t for t in tokens if not WELL_FORMED.match(t)})
    check("every braced token in the template is a well-formed placeholder",
          not bad, f"malformed: {bad}" if bad else f"{len(tokens)} tokens")

    # 2. the renderer's own pattern sees every token this audit sees
    seen = set(re.findall(R.PLACEHOLDER, tpl))
    unseen = sorted(set(tokens) - seen - set(bad))
    check("the renderer's pattern matches every well-formed token",
          not unseen,
          f"renderer blind to: {unseen}" if unseen else f"{len(seen)} tokens matched")

    # 3. every placeholder resolves to a key or a table
    keys = {t for t in tokens if not t.startswith("table:")}
    tbls = {t.split(":", 1)[1] for t in tokens if t.startswith("table:")}
    missing_k = sorted(keys - set(values))
    missing_t = sorted(tbls - set(tables))
    check("every value placeholder has a registry key", not missing_k,
          f"missing: {missing_k}" if missing_k else f"{len(keys)} keys")
    check("every table placeholder has a generator", not missing_t,
          f"missing: {missing_t}" if missing_t else f"{len(tbls)} tables")

    # 4. no registry key is dead, tested by perturbation rather than by
    #    pattern. A static scan cannot see a key built as v[f"within_{rule}"],
    #    so instead each key's value is changed and the manuscript re-rendered:
    #    if the output is byte-identical, nothing reads that key. This is a
    #    reachability test on the real render path, not a guess from source.
    baseline = R.render(tpl)

    def perturb(x):
        if isinstance(x, bool):
            return not x
        if isinstance(x, (int, float)):
            return type(x)(x + 7919)
        return str(x) + "_MUT7919"

    dead = []
    for k in sorted(values):
        mut = dict(values)
        mut[k] = perturb(values[k])
        try:
            out = R.render(tpl, values=mut, tables=MN.tables(mut))
        except Exception:
            continue          # a formatting error means the key IS read
        if out == baseline:
            dead.append(k)
    check("no registry key is unreachable from the manuscript", not dead,
          f"dead keys ({len(dead)}): {dead}" if dead
          else f"all {len(values)} keys reach the output")

    # 4b. the same rule for table generators. Checking only for unused KEYS
    #     let a stale "power_table" generator survive in the registry, still
    #     built from the 40-trial pilot after the template had moved to the
    #     400-trial surface: a second power table that no longer agreed with
    #     the printed one, invisible because nothing rendered it.
    dead_t = sorted(set(tables) - tbls)
    check("no table generator is unused by the manuscript", not dead_t,
          f"unused generators: {dead_t}" if dead_t
          else f"all {len(tables)} generators rendered")

    # 5. no placeholder resolves to something empty or placeholder-like
    empties = sorted(k for k in keys
                     if values[k] is None or str(values[k]).strip() in
                     ("", "nan", "None", "n/a", "NA", "-"))
    check("no placeholder resolves to an empty or missing value", not empties,
          f"empty: {empties}" if empties else "")
    empty_t = sorted(t for t in tbls if not str(tables[t]).strip())
    check("no table generator produces an empty table", not empty_t,
          f"empty: {empty_t}" if empty_t else "")

    # 6. nothing brace-delimited survives into the rendered file
    out = R.render(tpl)
    residue = ANY_BRACED.findall(out)
    check("rendered output contains no brace-delimited residue", not residue,
          f"residue: {residue[:5]}" if residue else "")

    # 7. the committed manuscript equals a fresh render
    check("MANUSCRIPT.md on disk equals a fresh render of the template",
          TARGET.exists() and TARGET.read_text() == out,
          "" if TARGET.exists() and TARGET.read_text() == out
          else "re-run `make manuscript`")
    return out


def mutate(label, text, expect):
    """Render mutated template; `expect` is 'reject' or 'accept'."""
    try:
        R.render(text)
        got = "accept"
    except Exception as e:
        got = "reject"
    ok = got == expect
    print(f"  {'PASS' if ok else 'FAIL'}  mutation: {label} — expected {expect}, got {got}")
    if not ok:
        FAILS.append(f"mutation: {label}")


def mutation_tests():
    tpl = TEMPLATE.read_text()
    anchor = "## 1. Introduction"
    assert anchor in tpl
    def inject(s):
        return tpl.replace(anchor, s + "\n\n" + anchor, 1)

    mutate("unknown lowercase key", inject("{{no_such_key_here}}"), "reject")
    mutate("unknown UPPERCASE key", inject("{{noSuchKeyHere}}"), "reject")
    mutate("unknown key, all caps", inject("{{NO_SUCH_KEY}}"), "reject")
    mutate("unknown table", inject("{{table:no_such_table}}"), "reject")
    mutate("single brace, not a placeholder", inject("{within_all}"), "accept")
    mutate("unclosed placeholder", inject("{{within_all"), "accept")
    mutate("stray closing braces", inject("some text}}"), "accept")
    # a real key must still render, in both letter cases that exist
    mutate("known lowercase key", inject("{{within_all}}"), "accept")
    mutate("known key with capitals", inject("{{kaiB_amp}}"), "accept")

    # an empty registry value must not render silently
    import manuscript_numbers as M
    real = M.collect
    try:
        M.collect = lambda: {**real(), "within_all": ""}
        R.PLACEHOLDER  # unchanged
        try:
            out = R.render(inject("{{within_all}}"))
            empty_ok = "FAIL"
        except Exception:
            empty_ok = "PASS"
    finally:
        M.collect = real
    # the renderer itself does not reject empties - check() 5 above does, on the
    # registry. Record which layer catches it rather than claiming both do.
    print(f"  NOTE  an empty registry value is caught by the registry check "
          f"(audit 5), not by the renderer: renderer {'rejects' if empty_ok=='PASS' else 'accepts'} it")


def main():
    print("== renderer audit")
    audit()
    print("\n== mutation tests")
    mutation_tests()
    print()
    if FAILS:
        print(f"RENDERER AUDIT FAILED: {len(FAILS)} check(s)")
        for f in FAILS:
            print(f"  - {f}")
        return 1
    print("RENDERER AUDIT OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
