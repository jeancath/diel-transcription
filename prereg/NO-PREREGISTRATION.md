# There is no pre-registration for this study

RULES.md §2 step 5 requires that `prereg/PREREGISTRATION.md` be frozen — hypotheses, models,
decision rules, SESOI, uncertainty estimator and sensitivity analyses — **before any confirmatory
model is fitted**, and never edited afterwards.

**That did not happen here.** This file exists instead of a pre-registration, and it is deliberately
not called one. No file in `prereg/` may be presented as a pre-registration of this work, and none
may be written later and dated earlier.

## What this means for the study's status

Under RULES.md §5, a decision rule adopted after the data are in view is **post hoc**. Several of
this study's central choices were made that way:

| Choice | When it was made |
|---|---|
| The marginal-preserving ortholog-identity null as the primary null | after the uniform-phase null was found to be inflated, which was after the first cross-lineage results |
| The within-strain pair as the calibration benchmark | after the cross-lineage effects were seen to be small |
| The rhythm-significance threshold for ortholog selection (q < 0.05, and the four-rule sensitivity around it) | during analysis; reported as a sensitivity sweep rather than a single pre-set rule |
| The power analysis and its jitter grid | during revision 3, to bound what the design could exclude |
| The functional-stratified null and the BRITE multi-assignment sensitivity | during revision 4-5 |
| Withdrawal of the §3.9 enrichment test | revision 6, after the test was found to hold a fitted alignment shift fixed across permutations |

Therefore: **the whole study is exploratory.** There is no confirmatory arm. The manuscript must not
describe any result as confirmatory, pre-specified or hypothesis-testing in the pre-registered sense,
and no p value in it carries pre-registered error control.

## What was fixed in advance, and what that is worth

Two things were genuinely fixed before the comparisons they govern, and both are verifiable in the
commit history rather than asserted here:

- The rhythm caller was validated on a *kaiABC*-null negative control and on *kaiBC* in wild type
  before any cross-lineage statistic was computed.
- The dataset registry — which series, which arms, which lineages — was committed before the
  cross-lineage comparison, so arms were not selected on their results.

Neither is a substitute for a pre-registration. They constrain two specific forms of selection, not
the analysis as a whole.

## What replaces it, imperfectly

What this study has instead is a frozen, mechanically enforced record of the decisions actually
taken:

- `results/deviations.md` and `DEVIATIONS.md` — numbered, dated changes across six revisions.
- `CLOSURE-REPORT.md` — three separate verdicts (content reproducibility, byte reproducibility,
  scientific validity), with what each does and does not establish.
- `results/ASSERTION_COVERAGE.csv` — every claim-level assertion, and whether a targeted mutation
  test proves it fires.
- `results/NUMBER_PROVENANCE.csv` — every statistic in the manuscript traced to the file that
  produced it.

These establish that the reported analysis is the analysis that ran. They do not establish that the
decision rules were chosen blind to the data, and nothing can establish that after the fact.

## For the next study

A pre-registration for the follow-up questions — the alignment-refitting null, the fine-stratum
power investigation, the Ferrari dataset reanalysis — is both possible and cheap, because this
study already fixes the estimator, the null and the power model those questions need. Write it on a
new branch before the data are touched (RULES.md §1: a new idea starts a new branch from `main`).
