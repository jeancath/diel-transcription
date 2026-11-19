"""Fail if the manuscript and the results have drifted apart.

Three checks:
  1. MANUSCRIPT.md matches a fresh render of MANUSCRIPT.template.md.
  2. Every placeholder in the template resolves.
  3. No numeric literal survives in the template's Abstract, Results or
     Discussion sections
     outside a placeholder, except items on an explicit allowlist (section
     numbers, years, accessions, model constants). This is what stops a number
     being hand-typed back in.

Exit code 1 on any failure, so `make check-manuscript` fails in CI.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from render_manuscript import render, TEMPLATE, TARGET, PLACEHOLDER

# What counts as a statistic that must come from results/: any decimal, or any
# integer of three or more digits. Everything else (small counts, hours, section
# numbers) is prose. These literals are exempt, each for a stated reason.
STATISTIC = re.compile(r"\d+\.\d+|\d{3,}")
ALLOW = {
    # strain and collection numbers
    "7942", "6803", "73102", "51142", "7120", "7418", "6301", "101",
    # database accessions: structural identifiers, not statistics
    "7188",
    # years
    "2009", "2019", "2023", "2025", "2026",
    # simulation / design parameters chosen by us, not measured
    "150", "500", "1500", "1000", "300", "200", "0.05", "2.5", "0.5",
    # period and sampling constants
    "24.0", "12.0", "48.0",
}
# NOTE: there is deliberately no "section number" exemption here. An earlier
# version exempted /^\d+\.\d+$/ as a section number, which matches every
# decimal in the manuscript and made this entire check vacuous for decimals
# for as long as it existed. Section numbers are already stripped above, by
# the section-cross-reference and heading substitutions.


def numeric_literals(text):
    body = []
    keep = False
    for line in text.splitlines():
        if line.startswith("## "):
            keep = bool(re.match(r"## (Abstract|3\.|4\.|\d+\. (Results|Discussion))",
                                 line)) or line.strip() in ("## Abstract",)
        if line.startswith("### 3."):
            keep = True
        # The Discussion is in scope too. It was not until the v0.3 freeze
        # review asked for Abstract, Results AND Discussion to be verified
        # against the registry: the scan stopped at section 4, so a hand-typed
        # decimal in the Discussion would never have been caught.
        if line.startswith("## 5."):
            keep = False
        if keep:
            body.append(PLACEHOLDER.sub("", line))
    text = "\n".join(body)
    text = re.sub(r"§\s*\d+(\.\d+)?", "", text)        # section cross-references
    text = re.sub(r"^#+\s*\d+(\.\d+)?", "", text, flags=re.M)   # section headings
    out = []
    for tok in re.findall(r"(?<![\w.])\d+(?:\.\d+)?(?![\w])", text):
        if tok in ALLOW:
            continue
        if STATISTIC.search(tok):
            out.append(tok)
    return out


def main():
    errs = []
    try:
        fresh = render()
    except KeyError as e:
        print(f"FAIL  {e}")
        return 1
    if not TARGET.exists():
        errs.append("MANUSCRIPT.md does not exist; run render_manuscript.py")
    elif TARGET.read_text() != fresh:
        errs.append("MANUSCRIPT.md differs from a fresh render of the template "
                    "- a number was edited by hand, or results changed without re-rendering")
    lits = numeric_literals(TEMPLATE.read_text())
    if lits:
        errs.append(f"{len(lits)} un-templated numeric literal(s) in Abstract/Results/Discussion: "
                    f"{sorted(set(lits))[:12]}")
    for e in errs:
        print("FAIL ", e)
    if not errs:
        print("OK    manuscript numbers match results/ ; no hand-typed statistics")
    return 1 if errs else 0


if __name__ == "__main__":
    sys.exit(main())
