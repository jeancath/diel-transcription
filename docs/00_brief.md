# Project brief

**Status: reconstructed, not contemporaneous.** RULES.md §2 step 3 asks for this file before any
analysis. It was not written then. What follows is a brief reconstructed on 19 Nov 2026 from the
committed record — `REPORT.md`, `LITERATURE-SEARCH.md`, `PRIOR-ART.md`, `DEVIATIONS.md` and the
commit history — and it is labelled as such so that nobody reads it as evidence that the design was
fixed in advance. Where the design changed while the data were in view, `prereg/NO-PREREGISTRATION.md`
says so and `results/deviations.md` records it.

## Question

Diel transcriptional programs in cyanobacteria were mapped principally in one strain,
*Synechococcus elongatus* PCC 7942, and are commonly discussed as if representative of the group.
Does the **identity and timing** of rhythmic orthologs transfer between cyanobacterial lineages?

Three claims are kept apart throughout, because the literature compresses them into one:

1. **Distribution of clock-protein homologues** — a statement about sequence presence.
2. **Functional conservation of the KaiABC oscillator** — a statement about a working,
   temperature-compensated, self-sustaining oscillator.
3. **Conservation of the transcriptional output program** — which genes peak, and when.

This study addresses (3), and only under light/dark driving. It cannot settle (2).

## Scope limit that defines the study

The data are diel: cells under light/dark cycles. A diel transcript rhythm can be driven by the
light/dark cycle, gated by an endogenous oscillator, or both, and these data cannot separate them.
Every conclusion is therefore about **diel transcriptional programs**, never about endogenous
circadian regulation. This is a construct-validity limit, not a wording preference.

## Hypotheses

| ID | Hypothesis | What would falsify it |
|---|---|---|
| H1 | Ortholog peak times agree between independent studies of the same strain | within-strain agreement no better than the phase-preserving null |
| H2 | Ortholog peak times agree between cyanobacterial lineages at a scale comparable to H1 | cross-lineage agreement far below the within-strain benchmark |
| H3 | A null expecting uniform phases inflates the apparent cross-species signal | its false-positive rate stays near nominal under realistic phase concentration |
| H4 | The time-point-permutation null used in the published Archaeplastida comparison is also inflated | its false-positive rate stays near nominal |
| H5 | Any residual cross-lineage correspondence survives conditioning on broad functional composition | the residual disappears once ortholog identity is permuted within functional classes |

H1–H4 were addressed as posed. H5 is answered only as a **sensitivity** result: the design cannot
distinguish loss of functional composition from a weak ortholog-specific signal from loss of power
in small strata.

## Data sources

7 public GEO series, 9 analysed arms, 6 lineages, 4 platforms. Sampling time recovered from
deposited sample metadata in every case. One *kaiABC*-null arm serves as a negative control.
The registry in `src/datasets.py` records, per arm, which series, which lineage and how sampling
time was recovered.

## Methods

- Cosinor rhythm detection with per-series intercepts, F-test, Benjamini-Hochberg FDR.
- Orthology by single-copy KEGG Orthology groups, and independently by Smith-Waterman reciprocal
  best hits (1175-1804 orthologs per pair).
- The central statistic: reduction in mean |ortholog phase difference| against a null, at each
  pair's best global time-origin shift.
- Five nulls, reported separately: uniform-phase, time-point permutation, marginal-preserving
  ortholog-identity permutation, operon-block, and functional-stratified.
- Power analysis by Monte Carlo over ortholog count, conserved fraction and phase jitter.

## Risks identified in the record

| Risk | How it was handled |
|---|---|
| The comparative design is not new — Ferrari et al. 2019 ran it for Archaeplastida and concluded the opposite | `PRIOR-ART.md`; the null-calibration analysis is framed as explaining the discrepancy, not as discovering the design |
| A null that ignores phase concentration manufactures significance | three nulls compared on simulated data with conservation fixed at zero |
| "Not detected" read as "absent" | enforced by a semantic assertion and a mutation test |
| A fitted alignment shift held fixed inside a permutation test | the affected test was withdrawn rather than reported with a caveat (§3.9) |
| Diel read as circadian | enforced by two mutation tests, one for the claim and one for the vocabulary |
| Numbers in the text drifting from the results files | every statistic rendered from `results/` and re-checked by `audit/number_check.py` |

## Phases as they actually ran

literature and prior-art sweep → data recovery and registry → rhythm caller and its negative
control → orthology → the within-strain benchmark → cross-lineage comparison → null calibration →
power analysis → functional stratification → six revisions → freeze at v0.3.

Note the ordering: the null-calibration and power work came **after** the first cross-lineage
results, in response to review. That is why this study is exploratory by default.
