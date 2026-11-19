# Deviation log, v0.3

What the closure sprint specified, what was actually done, and why. Entries
marked **deviation** depart from the brief; entries marked **correction** fix
something that was wrong before the sprint started.

---

## 1. Design-matched simulations

**Specified.** Complete the design-matched simulations, report them in full,
and do not attribute a simulated false-positive rate to the published analysis by name.

**Done as specified.** `results/n16_design_matched_nulls.csv`, 22 cells, 300
outer trials and 1000 inner permutations each, Wilson 95% intervals
throughout. Reported as two separate tables: (A) the primary
ortholog-identity null under the five sampling designs this manuscript
analyses, (B) a time-point-permutation null under three published sampling
designs. The banned phrase appears nowhere; §3.6 states explicitly what (B)
does not license, and the table generator carries the same statement as a
code comment so it survives a rewrite of the prose.

**Observation not anticipated by the brief.** Under the 7-point, 4 h design,
the inflation is not monotonic in phase concentration: 0.107 at 40%, 0.160 at
60%, 0.087 at 90%. Coarse sampling limits phase resolution, so extreme
concentration compresses the observed statistic as well as the null. This is
reported rather than smoothed over, and it weakens any attempt to read table
(B) as a correction factor.

## 2. Power surface

**Done as specified.** `results/n14_power_surface.csv` (45 cells) and
`results/n17_power_small_fractions.csv` (18 cells) together cover conserved
fraction 0, 10, 20, 40, 60, 80, 100% × jitter SD 1, 2, 4 h × n = 150, 500,
1500, at 400 outer trials and 400 inner permutations, each with a Wilson
interval. Size holds at every zero-conservation cell (0.030–0.065).

**Deviation.** The brief asked for one grid. It was run as two files because
the small-fraction cells were already complete when the sprint began;
re-running them inside one grid would have cost an hour and changed no number.
The manuscript table is assembled from both files, and the cell counts,
trial counts and seeds are identical.

## 3. Functional-stratified null

**Done as specified.** Stratum sizes, both BRITE levels, small-strata
sensitivity, exchangeability documented, pre-specified.
`results/n18_stratified_closure_A.csv` and `_B.csv`.

**Deviation from the expected conclusion.** The brief anticipated closing the
null. The result does not close it. At BRITE level B with classes smaller than
ten dropped, only 2 of 6 cross-lineage pairs remain significant — but that
filter removes 53–78 orthologs per pair and roughly halves two of them
(146 → 75, 171 → 93), so attenuation is confounded with loss of power. The
manuscript therefore claims robustness to *coarse* functional conditioning
only, and states that independence from functional composition is **not**
established at fine conditioning on this panel. A stronger claim was available
and is not made.

## 4. Semantic validation

**Done, and it found more than expected.** 27 assertions, all passing, that
test meaning rather than format: a clock-gene label must match its KEGG KO,
*kaiB* and *kaiC* must be different loci, arm and series counts are derived
from the registry rather than asserted, every accession must resolve to a
matrix that declares it, phases must lie in [0, 24), the *kaiABC*-null control
must stay at zero, and the manuscript must not contain claims stronger than
the results support.

Four assertions were added during the sprint in response to defects they then
caught (§5).

## 5. Defects found by the new checks

Each was found by automation, not by reading, and each is now covered by an
assertion that was confirmed to fail before the fix.

| defect | consequence had it shipped |
|---|---|
| **correction** — placeholder pattern was `[a-z0-9_:]+`, so any key with a capital letter was neither substituted nor reported missing | 21 lines of the manuscript carried literal `{{...}}`: the clock-gene control values, the leave-one-out gains, the whole stratified section and both its tables. `check_manuscript` reported OK throughout. |
| **correction** — `results/lit_ferrari_citing_screen.csv` fed the manuscript but no script produced it | the citation-screen counts were unverifiable, and `make clean` would have destroyed them |
| **correction** — `closure_report()` was defined after `if __name__ == "__main__"` | task 3's deliverable was unreachable from any make target |
| **correction** — `orthology` was not in `make all`, but `make clean` removes its output | the rebuild in task 5 would have failed at the manuscript step |
| **correction** — `diel-conservation/data` was a tracked symlink to an absolute scratch path | a dangling link in any clone |
| **correction** — "9–78 orthologs removed per pair" was typed by hand | the real range is 53–78 |
| **correction** — seven files in `results/` had no producer at all | exploration leftovers that nothing read; deleted |

## 6. Simulation reproducibility

**Not in the brief; added.** The four Monte Carlo jobs wrote their CSV only at
the end, and an interruption lost the whole run. Worse, each advanced one RNG
stream across cells, so a resumed run would have produced *different* numbers
from an uninterrupted one — silently. Each cell now checkpoints as it
completes and draws from a generator seeded from its own coordinates, so a
cell's value depends only on which cell it is.

## 7. Explicitly not done

The Ferrari reanalysis. The brief defers it to a separate project and no part
of it was begun: their expression data was not downloaded, and no statistic in
this manuscript is computed from it. Everything said about their procedure
comes from a design-matched simulation with truth fixed at zero, and the
manuscript says so wherever the subject arises.

---

## Revision 3: what the second clean-checkout rebuild exposed

The first clean rebuild failed on missing inputs. The second failed on
unreachable code, and chasing that failure found a defect in the check that was
supposed to guarantee the manuscript's numbers.

1. **`multiassign_report` was defined below its script's `__main__` block**, so
   `make strat` could not reach it. `n19_brite_multiassign.csv` existed in the
   working copy only because it had been produced interactively. This is the
   same defect shape as `closure_report` in the previous round: fixing the first
   instance by moving the entrypoint to the end of the file left no guard, and a
   later addition landed below it again. There is now an assertion that no
   top-level definition may follow a `__main__` block, with a mutation test.

2. **The no-hand-typed-statistics check was vacuous for every decimal.** It
   exempted tokens matching `^\d+\.\d+$` as "section numbers"; that pattern
   matches `41.9`, `0.13`, `2.9` and every other decimal in the manuscript. The
   check had therefore only ever caught bare three-or-more-digit integers, which
   is why `146` and `171` were caught in the previous round and nothing else
   was. With the exemption removed it immediately found 16 literals in the
   Abstract and Results.

3. **Three of those literal sets were stale, not merely un-templated.** The
   §3.4 power prose still quoted the 40-trial pilot (0.88-1.00, 0.53-0.78,
   0.07-0.65) while the Abstract quoted the 400-trial surface through the
   correct registry keys, so the two sections of the manuscript disagreed with
   each other. The time-point-permutation growth figures read 0.28->0.44 and
   0.76->0.88 where the current simulation gives 0.21->0.48 and 0.83->0.92.

4. **A stale second power table survived in the registry.** A `power_table`
   generator built from the pilot stayed after the template moved to the
   surface table; nothing rendered it, and the renderer audit checked for unused
   registry *keys* but not unused table *generators*. It now checks both.

5. **§3.9's control was computed on a different gene set from the core it was
   controlling** - and the conclusion reverses when that is fixed. The reported
   core requires a gene to be rhythmic in at least 4 of 5 lineages with a
   dispersion cutoff taken over all scorable genes; the control used "at least
   3" throughout. At the matched definition the core is 33 orthologs against
   20.4 expected under ortholog-identity shuffling (p < 0.002, 500 permutations,
   1.61x), where the mismatched comparison had given 48 against 41.9 (p = 0.13)
   and licensed the heading "No identifiable conserved gene set".

   The claim was **not** simply reversed. The null holds the fitted global
   alignment constant while the alignment was fitted on true ortholog pairings,
   so it is anti-conservative by an unquantified amount. §3.9 now reports a
   candidate set, states the alignment caveat as the binding limitation, and
   claims neither that a conserved core exists nor that the control refutes one.
   Known gap 13 records the null that would settle it, which is not run here.

6. **Four functions sat below their `__main__` blocks unreferenced.**
   `sd_null` was the §3.9 control and is now wired in and reported;
   `genome_similarity`, `decay` and `phase_test_offset_free` wrote nothing and
   no manuscript claim depended on them, so they were deleted rather than left
   looking like analyses that had been run. They remain in git history.

Assertion coverage is now reported rather than asserted: 26 assertion sites in
`validate_semantics.py` (35 run-time checks, since some loop over datasets), 11
mutation tests, 13 assertion sites proven to fire by a targeted mutation, 13
resting on review. `results/ASSERTION_COVERAGE.csv` carries the mapping. An
earlier round of this file claimed "29 assertions, each with a mutation test";
that was wrong on both halves and is corrected here.

Still explicitly not done: the Ferrari dataset reanalysis.

---

## Revision 4: the §3.9 test is withdrawn, not corrected

Revision 3 fixed a mismatched comparison in §3.9 and the corrected comparison came
out significant. Review rejected that as evidence, correctly, and the claim is
now withdrawn rather than reversed.

1. **§3.9 reports no p value and no fold enrichment.** The candidate set is
   selected through a global alignment shift fitted to the observed ortholog
   pairings, and that shift is held fixed across permutations. The null
   therefore does not reproduce the selection procedure that produced the
   observed set, so it cannot calibrate a statistic that depends on it: the
   test is anti-conservative by an amount this design cannot quantify. This is
   selection bias inside the test, not a conservative approximation, so a
   caveat beside the p value would not have been enough. The section now
   reports the set size, the dispersion cutoff and the membership criterion
   only, and states plainly that no valid test is available. The diagnostic
   counts remain in `results/n5_core_sd_null.csv` as a record.

   The registry keys for the null mean, the fold enrichment and both p values
   were **deleted**, not merely unused: the renderer audit fails on an unused
   key, so those numbers can no longer be reintroduced without the audit
   objecting. A semantic assertion fails if any p value or fold enrichment
   appears in §3.9, with a mutation test that puts `p < 0.002` back and
   requires that assertion to be the one that fails. The written disclaimer
   ("no p value and no fold enrichment") does not trip the check, because the
   check requires a digit beside the statistic - the trap an earlier
   disclaimer fell into.

   What would license a p value is recorded as gap 13: refit the alignment
   inside every permutation, or fix the shift by a rule independent of these
   data. Neither is run here, and the manuscript makes no claim that needs it.

2. **A provenance ledger for every printed number** (`src/number_provenance.py`,
   `results/NUMBER_PROVENANCE.csv`). For each of the 147 registry keys: the
   template lines it appears on, the result file it was computed from, the make
   target that produces that file, and a perturbation test proving the rendered
   manuscript changes when the key changes. The mapping is derived by static
   analysis of the registry source, not maintained by hand, because a
   hand-kept ledger is exactly what went stale in revision 3. All 147 keys trace
   to a result file and all 147 are proven to affect the output. Attribution
   is deliberately over-inclusive rather than under-inclusive: two keys list a
   second file that a shared local variable carried in.

   Writing the tracer exposed four of its own blind spots, each a real pattern
   in the registry: `ast.walk` is breadth-first, so statements inside `if`
   blocks were read after the top-level statements consuming them; keys bound
   to a loop variable; lists built with `.append()`; and names bound by a
   `for` target. Prose literals that are legitimately not statistics - section
   numbers, years, strain numbers, chosen simulation parameters - stay on the
   explicit allowlist in `check_manuscript.py`, each with a stated reason.
   The rule is provenance, not a ban on digits.

3. **Assertion coverage is reported, not claimed** (`src/coverage_report.py`,
   `results/ASSERTION_COVERAGE.csv`): 28 assertion sites, 35 run-time checks,
   12 mutation tests, 14 sites proven to fire by a targeted mutation, 14
   resting on review. The previous round's "29 assertions, each with a
   mutation test" was wrong on both halves. A mutation counts as coverage only
   when the validator fails **on the named assertion**, so a generic crash
   cannot be mistaken for it.

   One mutation test immediately found a weakness in itself rather than in the
   code: the BRITE declaration is made twice, in Methods and in Results, and
   removing one left the project compliant. The mutation now removes both.

4. **BRITE multi-assignment: "stable" is now defined.** §3.8 previously
   interpreted the *Prochlorococcus* pair in one paragraph and disowned its p
   value in the next. The interpretation is removed. "Stable under the
   assignment rule" is defined as: the count of significant cross-lineage
   pairs changes by at most one, and no pair's effect changes sign. It
   explicitly does not mean the identity of the significant pairs is stable
   (pairs do cross alpha) nor that effect sizes are stable (for small-effect
   pairs they are not). Methods and Results both state that first occurrence
   was fixed as the primary rule before any stratified p value was examined
   and that neither rule was chosen by its outcome, and the whole stratified
   comparison is labelled exploratory. Two assertions enforce those
   statements; one has a mutation test.

Still explicitly not done: the Ferrari dataset reanalysis.

---

## Revision 5: what the successful clean rebuild exposed

The clean-checkout rebuild at `680ecd5` succeeded: `MAKE_ALL_EXIT=0`, all 93
simulation cells complete (8 + 45 + 22 + 18), zero failed, zero skipped, no
`.failures.csv`. No source file changed during the run; `git status` in the
clean checkout shows modifications confined to `results/`, which the build
regenerates and the repository tracks.

Comparing that build against the pre-clean backup gave 26 files identical and
17 differences. Ten are new artifacts of features added after the backup was
taken (`PROVENANCE.json`, four simulation manifests, both `n18_*_last.csv`
sensitivity runs, `n19_brite_multiassign.csv`, `n3_rbh_params.json`,
`n5_core_sd_null.csv`). The remaining seven were worth the exercise:

1. **`n1_two_nulls.csv` lost two columns that no committed code has ever
   produced.** `effect_vs_marginal_pct` and `effect_vs_time_null_pct` appear in
   the committed CSV, but `git log -S` finds them in no version of any script:
   they came from an uncommitted interactive variant of
   `null_comparison.py`. Nothing reads them, so no claim is affected, and the
   clean build's file is the one the tracked code actually produces. This is a
   second instance of the pattern that produced `n19` interactively, found in
   data rather than in code.

2. **`make clean` deleted two files `make all` could not restore.**
   `n9_three_nulls.csv` is written by `three_nulls.py`, which is in no make
   target, and `determinism_report.json` comes from `tests-slow`. The pilot
   script and its output are deleted; `clean` now removes exactly what `all`
   regenerates, and the determinism artifacts have their own `clean-tests`.

3. **The registry had a dead fallback to that pilot.** If
   `n13_three_nulls_1000.csv` were ever missing, `manuscript_numbers.py` would
   silently publish 50-trial pilot numbers in place of the 1000-trial run. The
   branch could not fire in a clean build, because `make clean` had already
   deleted the pilot. It is removed and n13 is now required.

4. **"Every results file has a producing script" had never checked production.**
   It matched any `OUT / "name"` expression, which appears as often in a reader
   as in a writer: `n13_three_nulls_1000.csv` passed only because
   `manuscript_numbers.py` mentioned the path while *reading* it, and the
   assertion broke the moment that read was deleted. It now matches writing
   contexts only - `to_csv`, `Checkpoint`, `savefig`, `write_text` - and all 38
   result files pass against genuine writers.

5. **`n5_core_pathway_enrichment.csv` was not bit-reproducible.** Two pathways
   tied at p = 0.127155 swapped places between builds, because pandas'
   default sort is not stable. The sort is now `["p", "pathway"]` with
   `kind="mergesort"`.

6. **`n17_power_small_fractions.csv` differs numerically, as expected.** The
   committed file predates per-cell seeding; the differences are scattered in
   sign, mean +0.006, maximum 0.070 across 18 cells, against a standard error
   of about 0.032 for the difference of two independent estimates at 400
   trials. That is Monte Carlo resampling, not a change in the science. The
   quoted small-fraction power range moves from 0.11-0.58 to 0.08-0.60, and
   the manuscript takes the new values automatically because they are
   templated.

7. **`provenance.py` reported five outputs with no declared producer** - the
   three figures, `numbers.json` and `n19` - because its producer map lacked
   those prefixes. Declared.

Pipeline reproducibility and scientific validity remain separate verdicts, and
neither is claimed here: this round is reproducibility only.

Still explicitly not done: the Ferrari dataset reanalysis.

---

## Revision 6: the editorial conditions of the freeze

Review approved the v0.3 freeze subject to six editorial conditions. All six are
applied. The five that are statements in the manuscript are each enforced by an
assertion with a targeted mutation test, so a later edit that removes one fails
the build instead of passing quietly.

1. **§3.9 strictly exploratory.** The required sentence is in the text verbatim
   and asserted: the identified genes are an exploratory candidate set whose
   enrichment has not been established by a valid selection-adjusted null test.
   No p value or fold enrichment appears for that set, the registry keys for
   them stay deleted, and the candidate set is not mentioned in the Abstract.

2. **Functional stratification is sensitivity, not absence.** §3.8 now
   enumerates the three readings this panel cannot separate: the correspondence
   reflects functional composition; it is ortholog-specific but weak; or the
   test loses power under fine conditioning because small strata admit few
   distinguishable permutations. Loss of significance after stratification is
   stated to be evidence of sensitivity to functional conditioning, not
   evidence that ortholog-specific similarity is absent.

3. **Diel expression versus endogenous circadian regulation.** The Discussion
   states that this paper answers a question about diel transcriptional phase
   and cannot settle whether endogenous circadian regulation is conserved, with
   the reason: light-dark data confound oscillator gating with direct light
   response, and one lineage has adequate constant-light data. Two mutation
   tests guard it, one removing the sentence and one collapsing the vocabulary,
   because the distinction can be lost by word choice as easily as by claim.

4. **Review-based assertion sites are identified as such.** Coverage is now
   classified by consequence, not just counted: 31 assertion sites, 41 run-time
   checks, 22 mutation tests, 24 sites proven to fire by a targeted mutation, 7
   resting on documented review. Adding seven mutation tests during this round
   took the critical-but-untested count from 8 to 1. The remaining one, the RBH
   identity-case exclusion, resists data mutation because the exclusion lives
   in registry code rather than in any result file; it needs source mutation
   and is recorded as follow-up. `coverage_report.py` fails if any assertion
   site is unclassified, so a new critical check cannot default into a lenient
   bucket.

   Three of the seven new mutation tests found weaknesses in themselves rather
   than in the code, which is the argument for requiring a *named* assertion to
   fail: one removed a statement made twice and left the project compliant; one
   corrupted two fields and tripped a neighbouring assertion instead; one
   targeted a column holding a locus tag rather than a gene name.

5. **Content and byte reproducibility are separate verdicts.** Content passes.
   Byte-for-byte is incomplete by exactly one file: `rhythms_per_gene.csv.gz`,
   because gzip writes the current time into its header. The fix is one writer
   argument and is recorded as follow-up rather than bought with an hour-long
   rebuild for a header.

6. **Abstract, Results and Discussion verified against the registry.** This
   found a real gap: the hand-typed-statistics scan stopped at section 4, so a
   decimal in the Discussion would never have been caught. The Discussion is now
   in scope, verified by injecting `7.77%` and confirming the check fails.
   Extending it immediately flagged `7188`, from the ArrayExpress accession
   `E-MTAB-7188` - a structural identifier, now allowlisted with that reason.
   There is no separate Conclusion section; the Discussion closes the paper.

No new exploratory analysis was begun to satisfy these conditions. Deferred as
separate follow-up work, not as repairs to this paper: the alignment-refitting
null for §3.9, the fine-stratum power investigation, and the Ferrari dataset
reanalysis.
