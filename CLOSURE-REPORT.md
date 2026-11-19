# v0.3 Closure Report

Two verdicts, kept separate because they are different qualities and one does
not imply the other.

| | verdict |
|---|---|
| **Content reproducibility** | **passes.** Every number in the manuscript is regenerated from declared raw inputs by tracked source, and three independent full recomputations agree on every value. |
| **Byte-for-byte reproducibility** | **incomplete, by one file.** `rhythms_per_gene.csv.gz` is not byte-reproducible because gzip writes the current time into its header. Its content is identical, verified by reading through the compression. Every other output is byte-identical apart from the simulation manifests, which record their own wall time by design. |
| **Scientific validity** | **qualified. Not a blanket certification.** One section's test was withdrawn as invalid; eleven claims are narrowed; seven matters remain open and are listed at the end. |

The first two are deliberately separate rows. Identical content with a differing
file hash is a metadata fact; stating them as one verdict would either overclaim
byte reproducibility or understate content reproducibility. The fix for the
remaining case is one writer argument (`mtime=0`), recorded as follow-up rather
than bought with an hour-long rebuild for a header.

Commits. The authoritative build ran at **`327be5a`**. HEAD has since moved
(determinism manifests, and `compare_builds.py`, which `make all` does not
invoke). That the build's source equals HEAD's source is checkable, not a
claim to take on trust:

```bash
git diff --stat 327be5a HEAD -- diel-conservation/src/run_rhythms.py \
  diel-conservation/src/analyze.py diel-conservation/Makefile \
  diel-conservation/MANUSCRIPT.template.md
```

---

## 1. Rebuild — full log from a clean checkout

Three clean-checkout attempts. The first two **failed**, and are reported
because the failures are the evidence that this test was worth running.

| attempt | commit | result | what it exposed |
|---|---|---|---|
| 1 | `b0f1325` | **failed** at `nulls` | `fetch_data.py` never downloaded the 15 RBH inputs; `all` ran `nulls` before `orthology`; six undeclared files sat in the input cache |
| 2 | `680ecd5` | **failed** at `strat` | `multiassign_report` defined below its `__main__` block, so no make target could reach it; `n19` existed only because it had been produced interactively |
| 3 | `327be5a` | **`BUILD1_EXIT=0`** | - |

A fourth run was killed by a container restart mid-`nulls`. It was **not
resumed**: the clone's source was re-verified pristine and the build restarted
from `make clean`, so no output in the authoritative build predates the
restart.

Authoritative build: 46 declared inputs (26.7 MB, SHA-256 recorded), 49
outputs, every output with a declared producing target. `git status` in the
clean checkout shows no modification under `src/`, `Makefile` or the template
during the run.

`make clean` + `make all` is now a round trip. It was not: `n9_three_nulls.csv`
(from `three_nulls.py`, in no target) and `determinism_report.json` (from
`tests-slow`) were deleted by `clean` and never restored. The pilot is gone and
the determinism artifacts have their own `clean-tests`.

## 2. Determinism — uninterrupted versus interrupted and resumed

Re-run at the freeze commit, not carried over from an earlier round.

| | |
|---|---|
| cells compared | 8 |
| written before `SIGKILL` | 3 |
| recomputed after resume | 5 |
| checkpoint readable after `SIGKILL` | yes |
| resumed run equals uninterrupted run | **bit-identical, every cell and column** |
| negative control (shared RNG) | **differs in 25 values** |

The negative control is the point: it proves the comparison is capable of
failing. Under the original single-stream RNG design, a resumed run does not
reproduce an uninterrupted one, because skipped cells no longer consume draws.
Per-cell seeding from the cell coordinates fixes it.

Checkpoints are written to a temporary file, `fsync`ed and renamed, so an
interrupted write cannot leave a parseable truncated CSV. Each run writes a
manifest recording completed, failed and skipped cells against the expected
set, plus the seed, configuration, commit, dependency versions and wall time.

## 3. Simulations — every cell, with Monte Carlo intervals

| file | cells | failed | skipped | design |
|---|---|---|---|---|
| `n13_three_nulls_1000.csv` | 8/8 | 0 | 0 | 1000 outer x 200 inner |
| `n14_power_surface.csv` | 45/45 | 0 | 0 | 400 outer x 400 inner |
| `n16_design_matched_nulls.csv` | 22/22 | 0 | 0 | 300 outer x 1000 inner |
| `n17_power_small_fractions.csv` | 18/18 | 0 | 0 | 400 outer x 400 inner |

93 cells, no `.failures.csv`, all four manifests reporting
`completed == expected`. Every rate carries a Wilson 95% interval.

**The 93 cells were recomputed in full, three times, and agree bit-for-bit.**
`finish()` unlinks the partial checkpoint, so a re-run has nothing to resume
from: build 2 was a complete recomputation, not a skip. Wall times confirm it
(build 1: 1795.9 / 2714.0 / 1343.0 / 1095.1 s; build 2: 1804.5 / 2738.9 /
1334.7 / 1086.1 s). The third comparison is between the `680ecd5` and
`327be5a` clean builds, whose simulation source is byte-identical, run on
different container instances and different kernels (`v595`, `v653`).

**The coarse published design is not monotonic in phase concentration** and is
reported cell by cell, never as one inflation factor: 0.007 (0.002-0.024) at
uniform phases, 0.107 (0.077-0.147) at 40%, 0.160 (0.123-0.206) at 60%,
0.087 (0.060-0.124) at 90%. The 60% and 90% intervals do not overlap, so it is
not Monte Carlo noise. §3.6 gives the coarse-grid mechanism explicitly as an
explanation consistent with the numbers, **not** as one these runs establish,
because they do not record the observed and null locations separately within
that design.

These are simulations under synthetic phase distributions and study-like
sampling designs. They do not establish the behaviour of any published
analysis on its own data.

## 4. Manuscript — placeholder audit and source-to-claim mapping

| check | result |
|---|---|
| registry keys | 161 |
| keys traced to a result file | **161** |
| keys the rendered output provably depends on | **161** |
| as placeholders / as table-generator inputs | 150 / 11 |
| unresolved placeholders, any capitalisation | 0 |
| brace residue in the rendered output | 0 |
| keys resolving to empty | 0 |
| unused registry keys | 0 |
| unused table generators | 0 |
| hand-typed statistics in Abstract or Results | 0 |
| renderer mutation tests | 9, all firing as intended |

`results/NUMBER_PROVENANCE.csv` is the source-to-claim mapping: per key, the
template lines it appears on, the result file it was computed from, the
producing make target, and a perturbation test that re-renders the manuscript
with the key changed and requires the output to change. It is derived by static
analysis of the registry, not maintained by hand, because a hand-kept ledger is
what went stale. Attribution is deliberately over-inclusive: two keys list a
second file carried in by a shared local variable.

The audit re-derives the placeholder set independently of the renderer, because
a check that shares an implementation with the thing it checks proves nothing.
That is not hypothetical: the renderer's own pattern was lowercase-only, which
left 21 lines of the manuscript as literal `{{...}}` while `check_manuscript`
reported OK, since the same pattern both found and reported placeholders.

Structural numbers that are not statistics - section numbers, years, strain
numbers, database accessions, chosen simulation parameters - stay on an explicit
allowlist, each with a stated reason. The rule is provenance, not a ban on
digits. **Known limit:** the literal check treats any decimal, and any integer
of three or more digits, as a statistic; one- and two-digit integers in those
sections rest on the semantic assertions and on review (§6, gap 14).

**The Discussion was out of scope until this freeze.** The scan stopped at
section 4, so a hand-typed decimal there would never have been caught. Abstract,
Results and Discussion are now all scanned, verified by injecting `7.77%` into
the Discussion and confirming the check fails. Extending it immediately found a
real literal - `7188`, from the ArrayExpress accession `E-MTAB-7188` - which is
a structural identifier and is now allowlisted with that reason rather than
templated.

## 5. Semantics — assertion results and mutation tests

Reported as separate quantities, because they are separate quantities.

| | |
|---|---|
| assertion sites in `validate_semantics.py` | 31 |
| run-time checks (some loop over datasets) | 41 |
| failures in the authoritative build | **0** |
| mutation tests | 22 |
| assertion sites proven to fire by a targeted mutation | **24** |
| assertion sites resting on documented review | **7** |

**The 7 review-based sites are not experimentally validated.** They are
classified by what a silent failure of each would cost, with the action that
classification implies:

| category | mutation-tested | review only | action for the review-only sites |
|---|---|---|---|
| critical: biological identity | 7 | **0** | — |
| critical: statistical claim | 14 | **1** | add a targeted mutation test (follow-up) |
| structural / documentation | 3 | 5 | documented review may suffice |
| low impact: meta-check | 0 | 1 | review accepted |

The single critical site without a mutation test is *"reported RBH range
excludes the same-genome identity case"*. It resists data mutation because the
exclusion lives in the registry code rather than in any result file, so
corrupting an input cannot exercise it; testing it needs source mutation, which
is recorded as follow-up work. `results/ASSERTION_COVERAGE.csv` carries the
full mapping and fails the build if any assertion site is unclassified, so a
new critical check cannot land in a lenient bucket by default.

A mutation counts as coverage only when the validator fails **on the named
assertion**, so a generic crash cannot be mistaken for it. That rule earned its
keep three times during this closure: one mutation passed because the statement
it removed was made twice; one failed on a neighbouring assertion because it
corrupted two fields where one would do; and one targeted a column that holds a
locus tag rather than a gene name. All three were weaknesses in the tests, not
in the code, and were found only because the test demands a named failure.

An earlier round of this project claimed "29 assertions, each with a mutation
test". That was wrong on both halves.

Semantic provenance verified: gene identifier to gene name to KEGG KO; GEO
series to analysed arms (7 series, 9 arms); accession to retrieved record; RBH
pair to source sequences, with `n3_rbh_params.json` recording tool, matrix, gap
penalties, e-value, tie-breaking rule and the SHA-256 of all 21 proteomes.

**One assertion that had never worked.** "Every results file has a producing
script" matched any `OUT / "name"` expression, which appears as often in a
reader as in a writer. `n13_three_nulls_1000.csv` passed only because
`manuscript_numbers.py` mentioned the path while *reading* it, and the check
broke the moment that read was deleted. It now matches writing contexts only.

## 6. Reproducibility — first build against second, and against the backup

| comparison | identical | expected to differ | accounted for | **unexpected** |
|---|---|---|---|---|
| build 1 vs build 2 (same commit) | 47 | 0 | 0 | **0** |
| pre-clean backup vs build 1 | 27 | 2 | 19 | **0** |

Build-to-build runs in strict mode: nothing is allowed to differ except
volatile metadata. The 19 differences against the backup are each declared **by
name** in `BASELINE_DELTAS` with the change that caused them, and the tool
fails if a declared difference does not occur, so the table cannot rot in
either direction. Ten are new artifacts of validation work added after the
backup; the rest:

- **`n1_two_nulls.csv` lost two columns no committed code has ever produced.**
  `git log -S` finds `effect_vs_marginal_pct` in no version of any script, only
  inside committed CSVs: they came from an uncommitted interactive variant.
  Nothing reads them. The backup is therefore not a clean reference - it
  contains at least one file the tracked source cannot produce, which is the
  second instance of the pattern that produced `n19` interactively.
- **`n9_three_nulls.csv` removed**, with the dead registry fallback that would
  have silently published 50-trial pilot numbers in place of the 1000-trial run
  had `n13` ever been missing.
- **`n17_power_small_fractions.csv` differs numerically, legitimately.** The
  committed file predates per-cell seeding. 18 cells, differences scattered in
  sign, mean +0.006, maximum 0.070, against a standard error of about 0.032 for
  two independent estimates at 400 trials. The quoted small-fraction power range
  moves 0.11-0.58 to 0.08-0.60; the manuscript takes the new values because they
  are templated.
- **`n5_core_pathway_enrichment.csv` is now in the identical column.** It was
  not bit-reproducible: two pathways tied at p = 0.127155 swapped places
  between builds, because pandas' default sort is not stable.

`PROVENANCE.json` is compared on its meaning rather than its bytes: input
checksums must match exactly, and output checksums may differ only for a
declared reason. Result: inputs identical; of 49 outputs, the 4 simulation
manifests differ by recorded wall time, and `rhythms_per_gene.csv.gz` differs
because **gzip writes the current time into its header**. That file is not
byte-reproducible; its content is identical, verified by reading through the
compression. Declared rather than fixed, since the fix is a writer argument and
the content is verified.

## 7. Science — claims changed or narrowed after analysis

| was | is | why |
|---|---|---|
| "No identifiable conserved gene set" (§3.9) | a candidate set, **with no p value and no fold enrichment reported in either direction** | the control used a different membership criterion from the core it described; at the matched criterion it reverses, but that comparison is also invalid - the set is selected through an alignment shift fitted to the observed pairings which the null holds fixed, so it does not reproduce the selection procedure |
| *Prochlorococcus* "fails under every stratified variant" | **withdrawn**; no individual pair's stratified p value is interpreted | BRITE is not a partition; the assignment rule moves that pair from p = 0.129 to 0.009 at level A and 0.246 to 0.000 at level B |
| aggregate stratified result "stable" | **"stable" defined**: the count of significant pairs moves by at most one and no effect changes sign; explicitly **not** the identity of significant pairs, nor effect sizes | pairs do cross alpha under the alternative rule |
| power 0.88-1.00 / 0.53-0.78 / 0.07-0.65 (§3.4) | 0.87-1.00 / 0.53-0.74 / 0.08-0.60 | the prose quoted a 40-trial pilot while the Abstract quoted the 400-trial surface through the registry; the two sections disagreed |
| time-point null growth 0.28→0.44, 0.76→0.88 | 0.21→0.48, 0.83→0.92 | stale hand-typed literals |
| "9-78 orthologs per pair" | 53-78 | hand-typed and simply wrong |
| the published null "expects uniform phases" | a time-point-permutation null, substantially protected relative to a uniform null, and the appropriate object of criticism | v0.1 error |
| effect "absent" in *Prochlorococcus* pairs | "not detected, not shown to be zero" | p = 0.236 is not a claim of zero |
| one inflation factor for the published designs | each design and concentration read on its own row; no null model generally superior | the 7-point/4 h design is non-monotonic, with non-overlapping intervals |
| leave-one-out, phase-concentration correlations, functional stratification | labelled **exploratory**, over dependent pairs | ten pairwise comparisons share genomes and are not ten independent replications |
| priority | "to our knowledge", with the search's limits stated | one-level citation traversal, PubMed only, no preprint full text |

Enforced mechanically where wording can decay: 11 banned over-claim patterns, 5
required qualifications, and an assertion that §3.9 contains no p value or fold
enrichment - which looks for a digit beside the statistic, so the section's own
written disclaimer does not satisfy it. That trap is not hypothetical: an
earlier disclaimer tripped its own check.

---

## Editorial conditions of the freeze

Six conditions were set by review before v0.3 could be tagged. Each is applied,
and the five that are statements in the manuscript are enforced by an assertion
with a targeted mutation test, so that a later edit removing one fails the
build rather than passing silently.

| condition | how it is enforced |
|---|---|
| 1. §3.9 strictly exploratory; no p value, fold enrichment or confirmatory conserved-core claim | the required sentence is asserted verbatim; a second assertion rejects any p value or fold enrichment in §3.9, matching on a digit so the section's own disclaimer cannot satisfy it; the registry keys are deleted, so the unused-key rule blocks their return. Three mutation tests. |
| 2. functional stratification described as sensitivity to functional conditioning, not proof that ortholog-specific similarity is absent | §3.8 now enumerates the three readings this panel cannot separate - functional composition, weak ortholog-specific signal, and power loss under small strata - and states the sensitivity reading explicitly. Asserted; one mutation test. |
| 3. preserve the distinction between diel expression and endogenous circadian regulation | the Discussion states that this paper answers a diel-phase question and cannot settle whether endogenous circadian regulation is conserved. Asserted; two mutation tests, one removing the sentence and one collapsing the vocabulary. |
| 4. the review-based assertion sites identified as not experimentally validated | §5 above, with the risk classification and the per-category action. Reported, never implied to be validated. |
| 5. content-level and byte-level reproducibility kept distinct | two separate rows in the verdict table, and `PROVENANCE.json` is compared on its meaning: inputs must match exactly, outputs only for a declared reason. |
| 6. final Abstract, Results, Discussion and Conclusion verified against the frozen registry | the drift check requires the manuscript to equal a fresh render; the literal scan now covers the Discussion as well, proven by injection. There is no separate Conclusion section - the Discussion closes the paper. |

No new exploratory analysis was started to satisfy these conditions. Two
analyses are deliberately deferred and recorded as separate follow-up work: the
alignment-refitting null for §3.9, and the fine-stratum power investigation.

## What is not closed

1. **§3.9 has no valid test** (gap 13). Establishing enrichment needs the
   alignment refitted inside every permutation, or the shift fixed by a rule
   independent of these data. Not run; no claim depends on it.
2. **Fine-strata attenuation is confounded with power** - dropping small classes
   roughly halves two comparisons.
3. **One lineage has adequate constant-light data.** Everything else is diel
   under light-dark, so this is about diel transcriptional phase, not
   established endogenous circadian regulation.
4. **Seven assertion sites rest on documented review**, not on experiment
   (7 of 31). One is critical - the RBH identity-case exclusion - and resists
   data mutation because it lives in code rather than in a result file.
5. **Two- and one-digit integers in Abstract and Results** are not covered by
   the literal check (gap 14).
6. **`rhythms_per_gene.csv.gz` is not byte-reproducible** (gzip header).
7. **Only the simulations are resumable.** A container restart during
   `rhythms`, `analysis`, `nulls` or `orthology` requires a re-clean.

The Ferrari dataset reanalysis has not been started. It needs its own question,
hypotheses and protocol.
