"""Executable semantic assertions. The build fails on meaning, not just format.

The number-checker stops hand-typed statistics. It cannot stop a value that is
computed correctly from the wrong record - which is what produced every one of
the four errors found in v0.2: kaiB reported as kaiC, GEO series counted as
analysed arms, a Trichodesmium phase taken from the wrong gene row, and an RBH
range that silently included the within-strain identity case.

Each assertion below re-derives a fact from its source record and fails if the
pipeline disagrees. Run via `make validate` and as part of `make all`.
"""
import os
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from analyze import OUT
import orthology as O
from paths import GEO, KEGG
from datasets import DATASETS, BY_ID

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
RES = OUT                      # results/ , as every analysis script writes it

FAILS = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  — {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


def main():
    print("Semantic validation\n")

    # ---- 1. clock-gene identity: each row's KO must match its label
    KO_OF = {"kaiA": "K08480", "kaiB": "K08481", "kaiC": "K08482",
             "sasA": "K08479", "rpaA": "K10697", "cirA": "K29143"}
    cg = pd.read_csv(OUT / "clock_gene_control.csv")
    bad = cg[cg.apply(lambda r: KO_OF.get(r.clock_gene) != r.ko, axis=1)]
    check("clock-gene label matches its KEGG KO in every row", bad.empty,
          f"{len(bad)} mismatched rows" if len(bad) else f"{len(cg)} rows checked")

    # the specific confusion that occurred: kaiB and kaiC must be distinct rows
    wt = cg[cg.dataset == "SE7942_ito_WT"].set_index("clock_gene")
    for g in ("kaiB", "kaiC"):
        check(f"{g} present for the wild-type reference dataset", g in wt.index)
    if {"kaiB", "kaiC"} <= set(wt.index):
        check("kaiB and kaiC are different genes, not the same row repeated",
              wt.loc["kaiB", "gene"] != wt.loc["kaiC", "gene"],
              f"{wt.loc['kaiB','gene']} vs {wt.loc['kaiC','gene']}")

    # ---- 2. dataset bookkeeping: arms vs GEO series, derived from the registry
    ds = pd.read_csv(OUT / "dataset_summary.csv")
    check("every analysed arm appears in the dataset registry",
          set(ds.dataset) == {d.did for d in DATASETS},
          f"{len(ds)} arms")
    check("GEO series count is derived, not asserted",
          int(ds.gse.nunique()) == len({d.gse for d in DATASETS}),
          f"{ds.gse.nunique()} series / {len(ds)} arms")
    check("arms outnumber series only because of split arms",
          len(ds) >= ds.gse.nunique())

    # ---- 3. accessions resolve to files that exist and are non-empty
    missing = [d.gse for d in DATASETS
               if not (GEO / f"{d.gse}.txt.gz").exists()
               or (GEO / f"{d.gse}.txt.gz").stat().st_size == 0]
    check("every GEO accession resolves to a downloaded series matrix",
          not missing, f"missing: {sorted(set(missing))}" if missing else "all present")
    bad_acc = [d.gse for d in DATASETS if not re.fullmatch(r"GSE\d+", d.gse)]
    check("accession identifiers are well-formed", not bad_acc)

    # ---- 4. the series matrix actually belongs to the accession claimed
    import gzip
    wrong = []
    for d in {x.gse: x for x in DATASETS}.values():
        with gzip.open(GEO / f"{d.gse}.txt.gz", "rt", errors="replace") as fh:
            head = "".join(next(fh) for _ in range(40))
        if d.gse not in head:
            wrong.append(d.gse)
    check("each downloaded matrix declares the accession it was filed under",
          not wrong, f"mismatched: {wrong}" if wrong else "")

    # ---- 5. KO mapping sanity per dataset
    rh = pd.read_csv(OUT / "rhythms_per_gene.csv.gz")
    for did in ("SE7942_ito_WT", "SYN6803", "TERY", "MED4", "NPUN", "HALO_05"):
        sub = rh[rh.dataset == did]
        frac = (sub.ko.notna() & (sub.ko != "")).mean()
        check(f"{did}: a non-trivial fraction of genes carry a KO", frac > 0.10,
              f"{frac:.0%}")

    # ---- 6. RBH table excludes the within-strain identity case from its range
    rb = pd.read_csv(OUT / "n3_rbh_orthologs.csv")
    within = rb[rb.apply(lambda r: {r.a, r.b} == {"SE7942_ito_WT", "SE7942_vijayan"}, axis=1)]
    cross = rb[~rb.index.isin(within.index)]
    import manuscript_numbers as MN
    v = MN.collect()
    check("reported RBH range excludes the same-genome identity case",
          v["rbh_max"] == int(cross.n_rbh.max()) and v["rbh_min"] == int(cross.n_rbh.min()),
          f"reported {v['rbh_min']}–{v['rbh_max']}, cross-lineage "
          f"{int(cross.n_rbh.min())}–{int(cross.n_rbh.max())}, identity case "
          f"{int(within.n_rbh.iloc[0]) if len(within) else 'n/a'}")

    # ---- 7. phases are on the stated scale
    check("all fitted phases lie in [0, 24)",
          bool(((rh.phase_h >= 0) & (rh.phase_h < 24)).all()))
    check("all q values lie in [0, 1]",
          bool(((rh.q >= 0) & (rh.q <= 1)).all()))

    # ---- 8. the negative control is still negative
    kn = rh[(rh.dataset == "SE7942_ito_KAINULL")]
    check("kaiABC-null control yields no rhythmic genes at q<0.05",
          int((kn.q < 0.05).sum()) == 0, f"{int((kn.q < 0.05).sum())} of {len(kn)}")

    # ---- 9. manuscript claims are not stronger than the results
    man = (Path(__file__).resolve().parents[1] / "MANUSCRIPT.md").read_text()
    banned = {
        "proves": "no analysis here proves anything",
        "Ferrari et al. are wrong": "their data was never reanalysed",
        "demonstrates that clock rewiring": "rewiring is a hypothesis, not a result",
        "no conservation": "a weak conserved component is detected",
        # a simulated null's size under a design is not a property of anyone's
        # published result, and must not be named as though it were
        "ferrari's false-positive": "a design-matched simulation measures a "
                                    "procedure, not a published analysis",
        "ferrari's false positive": "a design-matched simulation measures a "
                                    "procedure, not a published analysis",
    }
    hits = [k for k in banned if k.lower() in man.lower()]
    check("manuscript avoids claims stronger than the results support",
          not hits, f"found: {hits}" if hits else "")
    check("manuscript distinguishes diel expression from circadian regulation",
          "diel" in man.lower() and "circadian regulation" in man.lower())
    check("manuscript states the Ferrari simulation is design-matched only",
          "design-matched" in man.lower() or "design-matched only" in man.lower()
          or "does not establish" in man.lower())

    # ---- claim-level wording audit
    #
    # Each entry is a specific over-claim this project is at risk of making,
    # not a style preference. The required phrases are the qualifications that
    # keep the corresponding claim inside what the analysis supports.
    low = man.lower()
    over = {
        "all cyanobacteria": "five lineages are not the phylum",
        "universal cyanobacterial": "no universality is established",
        "across cyanobacteria generally": "five lineages are not a general claim",
        "ten independent": "the ten pairs share genomes and are not independent",
        "independent replications": "the pairs are not replications",
        "first study to": "the literature search cannot establish priority",
        "no previous study exists": "absence of evidence from one search",
        "method is invalid": "their data was never reanalysed",
        "methodologically invalid": "their data was never reanalysed",
        "circadian output program is conserved": "these are diel time courses",
        "proves": "no analysis here proves anything",
    }
    hits = sorted(k for k in over if k in low)
    check("manuscript makes no claim wider than the design supports",
          not hits, f"found: {hits}" if hits else f"{len(over)} patterns checked")

    required = {
        "to our knowledge": "the priority claim must be hedged",
        "not independent": "pair dependence must be stated",
        "diel": "the LD datasets measure diel expression",
        "exploratory": "the pair-level and pile-up analyses are exploratory",
        "design-matched": "the published-design simulation must be labelled",
    }
    absent = sorted(k for k, why in required.items() if k not in low)
    check("manuscript carries the qualifications its claims require",
          not absent, f"missing: {absent}" if absent else
          f"{len(required)} qualifications present")

    # the leave-one-out and pile-up analyses must be labelled where they appear
    for term, label in (("leave-one", "leave-one-lineage-out"),
                        ("pile-up", "phase concentration"),
                        ("stratified null", "functional stratification")):
        if term in low:
            i = low.index(term)
            window = low[max(0, i - 1500):i + 1500]
            check(f"the {label} analysis is framed as exploratory or a sensitivity check",
                  any(w in window for w in
                      ("exploratory", "sensitivity", "descriptive", "cannot separate",
                       "do not choose", "not a formal")),
                  "")

    # the BRITE assignment rule is an analytical decision, and the manuscript
    # must say it was fixed before the p values were seen rather than chosen
    # for its outcome. A sensitivity analysis whose rule could have been picked
    # by result is not a sensitivity analysis.
    check("the BRITE assignment rule is declared as fixed in advance",
          "fixed as the primary rule before any stratified p value was examined" in low
          or "fixed as primary\nbefore any stratified p value was examined" in low
          or ("fixed as primary" in low and "before any stratified p value" in low),
          "")
    check("the manuscript states neither BRITE rule was chosen by its outcome",
          "not chosen for the result it gives" in low
          or "not selected by its outcome" in low, "")

    # ---- the rendered manuscript must contain no template residue
    #
    # This is checked by literal string search, not with the renderer's own
    # placeholder pattern. A bug in that pattern once left 21 lines of the
    # manuscript - including the clock-gene control values and both stratified
    # tables - as literal {{...}} text, substituted by nothing and reported as
    # missing by nothing, because the pattern did not match its own keys.
    residue = re.findall(r"\{\{[^}]{0,80}\}\}", man)
    check("rendered manuscript contains no unresolved template placeholder",
          not residue, f"found: {sorted(set(residue))[:5]}" if residue else "")

    # ---- every raw input the pipeline reads must be declared by the fetcher
    #
    # data/rbh/*.{faa,ft,gff}.gz were downloaded by hand and declared nowhere,
    # so `make all` from a clean checkout died on a missing file while every
    # check here passed. An input that only exists because someone fetched it
    # once is not a reproducible input.
    import fetch_data as FD
    declared = set()
    for g in FD.SERIES:
        declared.add(f"geo/{g}.txt.gz")
    for n in {**FD.SUPPL, **FD.GFF}:
        declared.add(f"geo/{n}")
    for pl in FD.PLATFORMS:
        declared.add(f"geo/{pl}.annot.txt")
    for key in FD.RBH_ASSEMBLIES:
        for suf in FD.RBH_SUFFIXES:
            declared.add(f"rbh/{key}.{suf}")
    present = set()
    for sub in ("geo", "rbh"):
        d = GEO.parent / sub
        if d.exists():
            for f in sorted(d.iterdir()):
                if f.is_file() and not f.name.endswith(".pkl"):
                    present.add(f"{sub}/{f.name}")
    undeclared = sorted(present - declared)
    check("every raw input present is declared by fetch_data.py",
          not undeclared, f"undeclared: {undeclared}" if undeclared
          else f"{len(present)} inputs, all declared")

    # ---- portability: no tracked file may point outside the repository
    #
    # diel-conservation/data was committed as a symlink into a session-specific
    # scratch directory. It worked here and would be a dangling link in any
    # clone, with paths.py then defaulting to a ./data that did not exist.
    import subprocess
    try:
        tracked = subprocess.run(["git", "ls-files", "-s"], cwd=ROOT.parent,
                                 capture_output=True, text=True, timeout=60).stdout
        links = []
        for line in tracked.splitlines():
            if line.startswith("120000"):
                rel = line.split("\t", 1)[1]
                tgt = (ROOT.parent / rel)
                if tgt.is_symlink():
                    dest = Path(os.readlink(tgt))
                    if dest.is_absolute() or ".." in dest.parts:
                        links.append(f"{rel} -> {dest}")
        check("no tracked file is a symlink outside the repository",
              not links, f"found: {links}" if links else "")
    except Exception as e:
        check("symlink portability check could run", False, str(e))

    # ---- provenance: every derived file must have a producer in src/
    #
    # A results file that no script writes cannot be regenerated, so a number
    # taken from it is unverifiable and `make clean` would destroy it. This is
    # not a formatting check: it fails on a file whose origin has been lost.
    # Only a WRITING context counts. An earlier version of this check matched
    # any `OUT / "name"` expression, which appears just as often in a reader:
    # n13_three_nulls_1000.csv passed because manuscript_numbers.py mentioned
    # the path while reading it, and the check broke the moment that read was
    # deleted. It had therefore never verified that anything produces a file.
    src_text = "\n".join(f.read_text() for f in sorted(SRC.glob("*.py")))
    WRITES = (
        r'to_csv\(\s*(?:OUT|FIG|RES)\s*/\s*f?"([^"]+)"',
        r'Checkpoint\(\s*(?:OUT|RES)\s*,\s*f?"([^"]+)"',
        r'savefig\(\s*(?:FIG|OUT)\s*/\s*f?"([^"]+)"',
        r'(?:write_text|write_bytes)\(\s*\)?\s*(?:OUT|FIG|RES)\s*/\s*f?"([^"]+)"',
        r'(?:OUT|FIG|RES)\s*/\s*f?"([^"]+)"\s*\)\s*\.write',
    )
    written = set()
    for pat in WRITES:
        for m in re.finditer(pat, src_text):
            written.add(m.group(1))
    pats = [re.compile("^" + re.sub(r"\{[^}]*\}", r"[^/]+",
                                    re.escape(w).replace("\\{", "{").replace("\\}", "}"))
                       + "$") for w in written]
    orphans = []
    for f in sorted(list(RES.glob("*.csv")) + list(RES.glob("*.csv.gz"))):
        if f.name.endswith(".partial.csv"):
            continue            # transient per-cell checkpoint, not a result
        if not any(pat.match(f.name) for pat in pats):
            orphans.append(f.name)
    check("every results file has a producing script in src/",
          not orphans, f"no producer: {orphans}" if orphans
          else f"{len(list(RES.glob('*.csv')))} csv files, all accounted for")

    # ---- reproducibility: `make all` must rebuild everything `make clean` removes
    mk = (ROOT / "Makefile").read_text()
    all_line = re.search(r"^all:((?:.*\\\n)*.*)$", mk, re.M).group(1)
    all_targets = set(all_line.replace("\\", " ").split())
    producers = {"fetch_data": "data", "run_rhythms": "rhythms", "analyze": "analysis",
                 "robustness": "robustness", "conserved_core": "nulls",
                 "three_nulls_fast": "simulations", "power_grid_full": "simulations",
                 "power_small_fractions": "simulations",
                 "design_matched_nulls": "simulations",
                 "stratified_null": "strat", "literature_screen": "literature",
                 "rbh": "orthology", "rbh_check": "orthology",
                 "figures": "figures",
                 "coverage_report": "validate",
                 "number_provenance": "validate"}
    reachable = set()
    for script, target in producers.items():
        if target in all_targets or target in ("data", "rhythms"):
            reachable.add(script)
    missing = sorted(t for s, t in producers.items()
                     if t not in all_targets and t not in ("data", "rhythms"))
    check("every producing target is reachable from `make all`",
          not missing, f"unreachable: {missing}" if missing else
          f"targets in all: {sorted(all_targets)}")

    # -- editorial conditions of the v0.3 freeze. Each of these three
    # sentences was required by review before the version could be tagged, and
    # each states a limit that a later edit could remove without changing any
    # number. They are checked here so that removing one fails the build.
    check("3.9 states the candidate set's enrichment is not established",
          "enrichment has not been established by a valid selection-adjusted "
          "null test" in " ".join(low.split()), "")
    check("stratification is framed as sensitivity, not absence of similarity",
          "not evidence that ortholog-specific similarity is absent"
          in " ".join(low.split()), "")
    check("the diel versus endogenous circadian scope limit is stated",
          "cannot settle the broader question of whether endogenous circadian "
          "regulation is conserved" in " ".join(low.split()), "")

    # -- the Abstract must carry the exploratory status itself. A reader who
    # sees only the Abstract must not take the study as confirmatory, so it is
    # not enough for prereg/ or the Discussion to say it: this looks inside the
    # Abstract alone.
    # no ^ anchors: a mutation test may hand the validator a flattened copy of
    # the manuscript, and the section must still be found in it
    m = re.search(r"## abstract\b(.*?)(?=## )", low, re.S)
    abstract = " ".join((m.group(1) if m else "").split())
    check("the Abstract states the study is exploratory, with no pre-registration",
          "this study is exploratory throughout" in abstract
          and "no pre-registration was frozen" in abstract,
          "" if abstract else "no Abstract section found")
    check("the Abstract states the diel versus endogenous circadian limit",
          "nothing here settles whether endogenous circadian regulation is "
          "conserved" in abstract, "")

    # -- scientific: 3.9 must report no inferential statistic.
    #
    # The identity-shuffle null for the candidate phase-conserved set holds a
    # global alignment fitted on the observed ortholog pairings fixed across
    # permutations, so it does not reproduce the selection procedure and cannot
    # calibrate the statistic. The decision is that no p value and no fold
    # enrichment for that set may appear in the manuscript at all, rather than
    # appearing beside a caveat. This looks for a NUMBER next to the statistic,
    # so the section's own written disclaimer ("no p value and no fold
    # enrichment") does not trip it - the same trap that a previous disclaimer
    # fell into. The pathway-enrichment q value is a different analysis and is
    # deliberately still allowed.
    man_f = ROOT / "MANUSCRIPT.md"
    if man_f.exists():
        mt = man_f.read_text()
        i = mt.find("### 3.9")
        sec = mt[i:mt.find("\n---", i)] if i >= 0 else ""
        offenders = []
        if re.search(r"\bp\s*[<=\u2264]\s*\d", sec):
            offenders.append("p value")
        if re.search(r"\d(?:\.\d+)?\s*(?:\u00d7|-fold|fold)", sec):
            offenders.append("fold enrichment")
        check("3.9 reports no p value or fold enrichment for the candidate set",
              bool(sec) and not offenders,
              f"inferential statistic present: {offenders}" if offenders
              else f"{len(sec)} chars of 3.9 checked")

    # -- semantic: the gene symbols named in 3.9 must actually be rows of the
    # core they describe. A prose gene list is a claim about specific source
    # rows, and no numeric check can see it: the numbers around it stay valid
    # if a symbol is wrong, stale or invented.
    core_f = OUT / "n5_conserved_core.csv"
    man = (ROOT / "MANUSCRIPT.md")
    if core_f.exists() and man.exists():
        ct = pd.read_csv(core_f, index_col=0)
        core = ct[(ct.n_rhythmic >= 4)
                  & (ct.phase_sd_h < ct.phase_sd_h.quantile(.25))]
        symbols = set()
        for nm in core.ko_name.fillna(""):
            head = nm.split(";")[0]
            for sym in head.split(","):
                sym = sym.strip()
                if sym:
                    symbols.add(sym.lower())
        m = re.search(r"candidate phase-conserved set(.{0,900})", man.read_text(), re.S)
        named = set(re.findall(r"\*([A-Za-z][A-Za-z0-9]{2,5})\*", m.group(1))) if m else set()
        absent = sorted(g for g in named if g.lower() not in symbols)
        check("every gene symbol named in 3.9 is a row of the core it describes",
              bool(named) and not absent,
              f"named but not in core: {absent}" if absent
              else f"{len(named)} symbols, all present among {len(symbols)}")

    # -- structural: nothing a make target calls may be defined after the
    # entrypoint. Python executes the __main__ block in file order, so a
    # function appended below it is unreachable from the command line even
    # though the module imports cleanly and every editor shows it defined.
    # Two separate defects of exactly this shape reached a clean rebuild.
    shadowed = []
    for f in sorted(SRC.glob("*.py")):
        lines = f.read_text().split("\n")
        try:
            i = next(k for k, l in enumerate(lines)
                     if l.startswith('if __name__ == "__main__":'))
        except StopIteration:
            continue
        for l in lines[i + 1:]:
            m = re.match(r"(?:def|class)\s+(\w+)", l)
            if m:
                shadowed.append(f"{f.name}:{m.group(1)}")
    check("no top-level definition follows the __main__ block",
          not shadowed, f"unreachable after entrypoint: {shadowed}" if shadowed
          else f"{len(list(SRC.glob('*.py')))} scripts checked")

    print()
    if FAILS:
        print(f"SEMANTIC VALIDATION FAILED: {len(FAILS)} assertion(s)")
        for f in FAILS:
            print(f"  - {f}")
        return 1
    print("SEMANTIC VALIDATION OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
