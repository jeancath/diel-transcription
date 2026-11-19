"""Prove the semantic assertions fail on deliberate errors.

A passing assertion means nothing unless it fires when the error it targets is
actually present. There are more assertions than mutation tests; which
assertions are covered is reported by src/coverage_report.py, not assumed. This script copies the project into a
temporary directory, introduces one specific error, runs the real validator
there as a subprocess, and requires both that it fails and that the named
assertion is the one that failed.

Mutations mirror errors this project actually made or nearly made: kaiB values
attributed to kaiC, a series counted as an arm, an accession that does not
match its matrix, an out-of-range phase, a manuscript claim stronger than the
result, a results file with no producer, and template residue in the output.

Run: python3 src/test_semantic_mutations.py
"""
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FAILS = []


def sandbox(tmp):
    """A working copy: real src and results, data symlinked, manuscript copied."""
    d = Path(tmp) / "diel-conservation"
    (d / "results").mkdir(parents=True)
    shutil.copytree(ROOT / "src", d / "src",
                    ignore=shutil.ignore_patterns("__pycache__"))
    for f in list((ROOT / "results").glob("*.csv")) + \
             list((ROOT / "results").glob("*.csv.gz")) + \
             list((ROOT / "results").glob("*.json")):
        shutil.copy2(f, d / "results" / f.name)
    (d / "figures").mkdir(exist_ok=True)
    for n in ("MANUSCRIPT.md", "MANUSCRIPT.template.md", "Makefile"):
        if (ROOT / n).exists():
            shutil.copy2(ROOT / n, d / n)
    (d / "data").symlink_to((ROOT / "data").resolve())
    # a git dir so the symlink/commit checks can run; the parent is enough
    return d


def run_validator(d):
    r = subprocess.run([sys.executable, "src/validate_semantics.py"], cwd=d,
                       capture_output=True, text=True, timeout=900)
    return r.returncode, r.stdout + r.stderr


def case(label, mutate, expect_assertion):
    """Mutate a sandbox, run the validator, require the named assertion to fail."""
    with tempfile.TemporaryDirectory() as tmp:
        d = sandbox(tmp)
        try:
            mutate(d)
        except Exception as e:
            print(f"  FAIL  {label} — could not apply mutation: {e}")
            FAILS.append(label)
            return
        rc, out = run_validator(d)
        failed_lines = [l.strip() for l in out.splitlines() if l.strip().startswith("FAIL")]
        hit = any(expect_assertion.lower() in l.lower() for l in failed_lines)
        if rc != 0 and hit:
            print(f"  PASS  {label}")
        elif rc == 0:
            print(f"  FAIL  {label} — validator PASSED a corrupted project")
            FAILS.append(label)
        else:
            print(f"  FAIL  {label} — failed, but not on '{expect_assertion}'")
            print(f"        actual: {failed_lines[:3]}")
            FAILS.append(label)


# ---- mutations

def swap_kaib_kaic(d):
    f = d / "results" / "clock_gene_control.csv"
    t = pd.read_csv(f)
    b = t.clock_gene == "kaiB"
    c = t.clock_gene == "kaiC"
    t.loc[b, "clock_gene"] = "kaiC_tmp"
    t.loc[c, "clock_gene"] = "kaiB"
    t.loc[t.clock_gene == "kaiC_tmp", "clock_gene"] = "kaiC"
    t.to_csv(f, index=False)


def break_accession(d):
    f = d / "src" / "datasets.py"
    s = f.read_text()
    m = re.search(r'"(GSE\d+)"', s)
    s = s.replace(m.group(1), "GSE99999999", 1)
    f.write_text(s)


def phase_out_of_range(d):
    f = d / "results" / "rhythms_per_gene.csv.gz"
    t = pd.read_csv(f)
    t.loc[t.index[0], "phase_h"] = 25.7
    t.to_csv(f, index=False, compression="gzip")


def q_out_of_range(d):
    f = d / "results" / "rhythms_per_gene.csv.gz"
    t = pd.read_csv(f)
    t.loc[t.index[0], "q"] = 1.5
    t.to_csv(f, index=False, compression="gzip")


def kainull_rhythmic(d):
    f = d / "results" / "rhythms_per_gene.csv.gz"
    t = pd.read_csv(f)
    m = t.dataset == "SE7942_ito_KAINULL"
    idx = t[m].index[:5]
    t.loc[idx, "q"] = 0.001
    t.to_csv(f, index=False, compression="gzip")


def overclaim(d):
    f = d / "MANUSCRIPT.md"
    f.write_text(f.read_text() + "\n\nThis analysis proves the clock output is conserved.\n")


def orphan_result(d):
    (d / "results" / "zz_untraceable_numbers.csv").write_text("a,b\n1,2\n")


def template_residue(d):
    f = d / "MANUSCRIPT.md"
    f.write_text(f.read_text() + "\n\nEffect size was {{within_all}}%.\n")


def drop_circadian_regulation_term(d):
    """Remove every mention of circadian regulation as a distinct concept.

    The assertion requires the manuscript to keep diel expression and circadian
    regulation apart; collapsing the vocabulary is how that distinction is lost
    in practice, not by an explicit claim.
    """
    f = d / "MANUSCRIPT.md"
    t = f.read_text()
    assert "circadian regulation" in t.lower(), "no 'circadian regulation' to remove"
    import re as _re
    f.write_text(_re.sub("circadian regulation", "rhythmic output", t,
                         flags=_re.I))


def drop_design_matched_label(d):
    """Strip the design-matched label from the published-design simulation."""
    f = d / "MANUSCRIPT.md"
    t = f.read_text()
    assert "design-matched" in t.lower(), "no design-matched label to remove"
    import re as _re
    t = _re.sub("design-matched", "equivalent", t, flags=_re.I)
    t = t.replace("does not establish", "shows")
    f.write_text(t)


def insert_wide_claim(d):
    """Insert a claim wider than the design supports."""
    f = d / "MANUSCRIPT.md"
    f.write_text(f.read_text() +
                 "\n\nThis proves the circadian output program is conserved "
                 "across all cyanobacteria.\n")


def drop_required_qualification(d):
    """Remove the pair-dependence qualification the claims require.

    It is stated more than once, so every occurrence has to go: this is the
    second mutation to have found a weakness in itself rather than in the code.
    """
    f = d / "MANUSCRIPT.md"
    t = f.read_text()
    flat = " ".join(t.split())
    assert "not independent" in flat, "pair-dependence qualification not found"
    f.write_text(flat.replace("not independent", "distinct"))


def drop_brite_outcome_statement(d):
    """Remove only the statement that neither BRITE rule was chosen by outcome."""
    f = d / "MANUSCRIPT.md"
    t = f.read_text()
    n = 0
    for old in ("not chosen for the result it gives",
                "Neither rule was selected by its outcome"):
        if old in t:
            t = t.replace(old, "used throughout")
            n += 1
    assert n >= 1, "BRITE outcome statement not found"
    f.write_text(t)


def identical_kaib_kaic(d):
    """Give kaiB and kaiC identical values, as if one row were reused for both."""
    import pandas as _pd
    f = d / "results" / "clock_gene_control.csv"
    t = _pd.read_csv(f)
    # `gene` holds a locus tag; the gene NAME is in `clock_gene`
    col = "clock_gene"
    b = t[t[col].astype(str).str.fullmatch("kaiB", case=False, na=False)]
    c = t[t[col].astype(str).str.fullmatch("kaiC", case=False, na=False)]
    assert len(b) and len(c), f"kaiB/kaiC rows not found in {col}"
    # Only the locus tag is copied. Copying `ko` as well would break the
    # label-to-KO assertion instead, and a mutation has to fail the assertion
    # it targets, not merely fail somewhere.
    ref = t[(t.clock_gene.str.fullmatch("kaiB", case=False, na=False))
            & (t.dataset == "SE7942_ito_WT")]
    tgt = t[(t.clock_gene.str.fullmatch("kaiC", case=False, na=False))
            & (t.dataset == "SE7942_ito_WT")]
    assert len(ref) and len(tgt), "wild-type kaiB/kaiC rows not found"
    t.loc[tgt.index, "gene"] = t.loc[ref.index[0], "gene"]
    t.to_csv(f, index=False)


def drop_an_arm(d):
    """Delete one analysed arm from the dataset registry.

    The arm count, the series count and every per-arm number in the manuscript
    are derived from this file, so losing a row is how a dataset would silently
    disappear from the panel.
    """
    import pandas as _pd
    f = d / "results" / "dataset_summary.csv"
    t = _pd.read_csv(f)
    assert len(t) > 1, "nothing to drop"
    t.iloc[1:].to_csv(f, index=False)


def _drop(d, needle, label):
    f = d / "MANUSCRIPT.md"
    t = f.read_text()
    flat = " ".join(t.split())
    assert needle in flat, f"{label} not found in the rendered manuscript"
    # the sentence is line-wrapped, so remove it on the flattened text and
    # write the flattened result back; the validator flattens before matching
    f.write_text(flat.replace(needle, "", 1))


def drop_39_exploratory(d):
    """Remove the sentence review required in 3.9 before v0.3 could be tagged."""
    _drop(d, "enrichment has not been established by a valid selection-adjusted "
             "null test", "3.9 exploratory sentence")


def drop_strat_sensitivity(d):
    """Remove the sentence forbidding 'absence of similarity' readings of 3.8."""
    _drop(d, "not evidence that ortholog-specific similarity is absent",
          "stratification sensitivity sentence")


def drop_circadian_scope(d):
    """Remove the diel versus endogenous circadian scope limit."""
    _drop(d, "cannot settle the broader question of whether endogenous circadian "
             "regulation is conserved", "circadian scope sentence")


def drop_abstract_exploratory(d):
    """Remove the Abstract's statement that the study is exploratory.

    A reader who sees only the Abstract must not take the study as
    confirmatory. The statement exists in prereg/ and in the Discussion too,
    so the assertion is scoped to the Abstract alone; this mutation proves that
    scoping works.
    """
    _drop(d, "This study is exploratory throughout. No pre-registration was frozen",
          "Abstract exploratory statement")


def drop_abstract_circadian_limit(d):
    """Remove the Abstract's diel versus endogenous circadian limit."""
    _drop(d, "nothing here settles whether endogenous circadian regulation is "
             "conserved", "Abstract circadian limit")


def unfixed_brite_rule(d):
    """Delete the declaration that the BRITE assignment rule was fixed in advance.

    A sensitivity analysis whose assignment rule could have been chosen after
    seeing the p values is not a sensitivity analysis, and the per-pair results
    here do move across alpha when the rule is swapped.
    """
    f = d / "MANUSCRIPT.md"
    t = f.read_text()
    # the declaration is made twice, in Methods and in Results; a mutation that
    # removes only one leaves the project compliant, so both must go
    n = 0
    for old in ("fixed as the primary rule before any stratified p value\nwas examined",
                "fixed as primary\nbefore any stratified p value was examined",
                "not chosen for the result it gives",
                "Neither rule was selected by its outcome"):
        if old in t:
            t = t.replace(old, "used throughout")
            n += 1
    assert n >= 2, f"expected the BRITE declarations in the manuscript, found {n}"
    f.write_text(t)


def pvalue_in_39(d):
    """Put the withdrawn p value back into 3.9.

    The identity-shuffle null there is anti-conservative by an unquantified
    amount, so no p value for that set may be published at all. This is the
    mutation that proves the rule is enforced and not merely written down.
    """
    f = d / "MANUSCRIPT.md"
    t = f.read_text()
    old = "**We therefore report no p value and no fold enrichment for\nthis set.**"
    assert old in t, "3.9 disclaimer not found in the rendered manuscript"
    f.write_text(t.replace(old, "The excess is significant (p < 0.002, 1.61-fold).", 1))


def def_after_main(d):
    """Append a function below a script's entrypoint, where nothing can call it.

    Two real defects had this exact shape: closure_report and multiassign_report
    were both defined after their __main__ block, so `make strat` could not
    reach them, while importing the module and calling them by hand worked.
    """
    f = d / "src" / "stratified_null.py"
    f.write_text(f.read_text() + "\n\ndef unreachable_helper():\n    return 1\n")


def wrong_core_gene(d):
    """Name a gene in 3.9 that is not a row of the core it describes."""
    f = d / "MANUSCRIPT.md"
    t = f.read_text()
    assert "*fabI*" in t, "expected *fabI* in the rendered manuscript"
    f.write_text(t.replace("*fabI*", "*nifH*", 1))


def main():
    print("== semantic mutation tests")
    case("kaiB and kaiC labels swapped against their KOs",
         swap_kaib_kaic, "clock-gene label")
    case("a GEO accession that resolves to nothing",
         break_accession, "accession")
    case("a fitted phase outside [0, 24)",
         phase_out_of_range, "phases lie in")
    case("a q value outside [0, 1]",
         q_out_of_range, "q values lie in")
    case("the kaiABC-null control reporting rhythmic genes",
         kainull_rhythmic, "kaiABC-null")
    case("a manuscript claim stronger than the results support",
         overclaim, "claims stronger")
    case("a results file no script produces",
         orphan_result, "producing script")
    case("template residue in the rendered manuscript",
         template_residue, "unresolved template placeholder")
    case("a function defined below its script's entrypoint",
         def_after_main, "follows the __main__ block")
    case("a gene symbol in 3.9 that is not in the core",
         wrong_core_gene, "named in 3.9")
    case("a p value reinstated in 3.9",
         pvalue_in_39, "no p value or fold enrichment")
    case("the BRITE rule no longer declared as fixed in advance",
         unfixed_brite_rule, "fixed in advance")
    case("3.9's exploratory-only sentence removed",
         drop_39_exploratory, "enrichment is not established")
    case("3.8's sensitivity-not-absence sentence removed",
         drop_strat_sensitivity, "sensitivity, not absence of similarity")
    case("the diel/circadian scope limit removed",
         drop_circadian_scope, "endogenous circadian scope limit")
    case("the Abstract's exploratory statement removed",
         drop_abstract_exploratory, "Abstract states the study is exploratory")
    case("the Abstract's diel/circadian limit removed",
         drop_abstract_circadian_limit, "Abstract states the diel versus")
    case("diel and circadian regulation collapsed into one term",
         drop_circadian_regulation_term, "distinguishes diel expression")
    case("the design-matched label stripped",
         drop_design_matched_label, "design-matched only")
    case("a claim wider than the design inserted",
         insert_wide_claim, "no claim wider than")
    case("a required qualification removed",
         drop_required_qualification, "qualifications its claims require")
    case("the BRITE outcome-neutrality statement removed",
         drop_brite_outcome_statement, "chosen by its outcome")
    case("kaiB and kaiC given identical values",
         identical_kaib_kaic, "different genes")
    case("an analysed arm dropped from the dataset registry",
         drop_an_arm, "every analysed arm appears")
    print()
    if FAILS:
        print(f"MUTATION TESTS FAILED: {len(FAILS)} assertion(s) did not catch "
              f"their error")
        for f in FAILS:
            print(f"  - {f}")
        return 1
    print("MUTATION TESTS OK: every assertion fired on the error it targets")
    return 0


if __name__ == "__main__":
    sys.exit(main())
