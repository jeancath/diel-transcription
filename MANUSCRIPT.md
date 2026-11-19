# Within-strain reproducibility sets the scale for cross-species comparison of diel transcriptional programs

**DRAFT v0.3 — not submitted, not peer reviewed, not for preprint deposition.**
Revised in response to revision of v0.1 and v0.2; changes are listed in §7.
Every statistic is generated from `results/` by `src/render_manuscript.py`;
`make check-manuscript && make validate` must pass before this file is read as current.
Code and data: https://github.com/jeancath/diel-transcription (public release). The working
repository is private; the public release is an allowlisted snapshot of the same source and
results.

---

## Abstract

Diel transcriptional programs in cyanobacteria were characterised principally in one strain,
*Synechococcus elongatus* PCC 7942, and are commonly discussed as though representative of the
group. Whether the identity and timing of rhythmic orthologs transfers between cyanobacterial
lineages has not been measured. We reanalysed 7 public diel and circadian transcriptome GEO series
(9 analysed arms) from 6 lineages and four platforms, recovering sampling time from deposited sample
metadata in every case. Rhythm detection was validated on a *kaiABC*-null negative control (0
rhythmic genes of 2490 at FDR < 0.05) and on *kaiBC* in wild type (3.07 log₂ peak-to-trough,
dusk-phased at 11.3/11.4 h, q = 1.3e-09/8.8e-09).

We calibrate cross-lineage comparison against an internal benchmark that cross-species studies have
lacked: two independent studies of the same strain, from different laboratories on different
platforms. Within that strain, ortholog peak times agree to a median of 1.67 h, and the mean
|phase difference| falls 4.15 h (75.0%) below a null preserving each dataset's phase distribution.
Across lineages the same statistic falls 0.07–0.55 h (median 5.0%, maximum 12.3%). The gap holds under
four ortholog-selection rules including none, under two orthology assignments (single-copy KEGG
Orthology; Smith–Waterman reciprocal best hits, 1175–1804 orthologs per pair), and under nulls
preserving operon structure.

Power analysis across sample sizes (150–1500 orthologs), conservation levels and jitter widths shows the
design detects conservation of 40% of orthologs with power 0.87–1.00 under tight phase jitter
(1–2 h) but only 0.53–0.74 under loose jitter (4 h), and reaches only 0.08–0.60 at 10–20%. Under tight-jitter assumptions we therefore exclude timing conservation involving 40% of
orthologs; that exclusion is substantially weaker under 4 h jitter, and we cannot resolve the
residual. We do not claim to exclude a majority of the program: 40% is not a majority, and the
bound is conditional on the simulation's assumptions. The residual is also unevenly distributed: it is carried
disproportionately by pairs involving *S. elongatus* (median 10.1% versus 1.2% for pairs without it)
and is not detected in pairs involving *Prochlorococcus* — not detected, not shown to be zero —
and effect size correlates negatively, in an exploratory analysis over dependent pairs, with how
concentrated the marginal phase distributions are (ρ = -0.62 to -0.55), so sensitivity and biology
are not separable here.

Conditioning additionally on broad functional composition — permuting ortholog identity only within
KEGG BRITE classes — leaves the within-strain benchmark at 73.5% and keeps
5 of 6 cross-lineage pairs significant at the coarse class level, but only
2 of 6 under fine classes with small strata removed, where the
permutable set also shrinks substantially. Individual pairs' stratified p values are not
interpretable: 38.6% of KEGG Orthology groups carry more than one BRITE class, and
resolving that differently moves 2 of 6 pairs across α = 0.05 without
changing the aggregate. Phase correspondence alone is therefore not sufficient evidence of a
conserved regulatory program.

Separately, we compare three nulls on simulated expression with conservation fixed at zero. A null
expecting uniform phases has a false-positive rate of 0.80–1.00 under realistic phase concentration;
a time-point-permutation null of the kind used in a published Archaeplastida comparison is partially
but not fully protected (0.21–0.92); an ortholog-identity-permutation null preserving the marginals
holds 0.04–0.06. Every dataset examined here lies in the regime where the first two are inflated.

This study is exploratory throughout. No pre-registration was frozen, and several central analysis
choices — the primary null, the within-strain calibration benchmark, the ortholog-selection
threshold, the power model and the functional stratification — were made with the data in view, so
no p value reported here carries pre-specified error control. The data are also diel: cells under
light/dark cycles, in which a transcript rhythm can be driven by the light/dark cycle, gated by an
endogenous oscillator, or both. These data cannot separate those, and nothing here settles whether
endogenous circadian regulation is conserved.

---

## 1. Introduction

Three claims about the cyanobacterial clock are routinely compressed into one and should be kept
apart:

1. **Distribution of clock-protein homologues.** KaiB and KaiC homologues occur well beyond
   cyanobacteria [Schmelling et al. 2017]. This is a statement about sequence presence.
2. **Functional conservation of the KaiABC oscillator.** That a genome carries homologues does not
   establish a functioning, temperature-compensated, self-sustaining oscillator; demonstrations of
   function remain concentrated in a few organisms.
3. **Conservation of the downstream transcriptional program.** Whether the same orthologs cycle, at
   the same time of day, in different lineages.

Claim 1 is well supported. Claim 3 is the subject of this paper and, as far as our search
established, has not been tested within cyanobacteria.

A further distinction governs what our data can support. **Circadian regulation** is endogenous
timekeeping that persists in constant conditions; **diel expression** is variation across a
light/dark cycle; **light-responsive expression** is direct response to illumination. Free-running
(constant-light) data separate the first from the third; light/dark data do not. Six of our eight
time courses are light/dark. Our primary claims are therefore about **diel transcriptional
programs**, and we treat clock rewiring as one hypothesis among several for the differences we
measure, not as a demonstrated cause.

Two prior studies frame the work. Ferrari et al. [2019] compared diurnal transcriptomes across nine
members of Archaeplastida — a clade comprising glaucophytes, red and green algae and land plants,
with *Synechocystis* sp. PCC 6803 as an outgroup — using one-to-one orthologs, mean |phase
difference| and a global shift search, and concluded that diurnal programs are "remarkably similar"
despite large phylogenetic distance. They did not compare cyanobacteria with each other. Brooks et
al. [2023] analysed 57 mouse-liver time series and found core-clock phase highly consistent across
studies while rhythmic gene-set membership was not; they attributed the discrepancy to a combination
of technical library-preparation differences, biological and design heterogeneity, and the
statistics of per-study thresholding, and showed that joint cross-study modelling recovers
consistency that per-study analysis misses.

The implication we take from Brooks et al. is that a conservation measurement needs a reference
scale for what "the same program, measured twice" yields under the identical statistic. We supply
one and interpret cross-lineage conservation relative to it.

---

## 2. Methods

### 2.1 Datasets

7 GEO series, 9 analysed arms, from 6 lineages (Table 1), 2009–2025, across three microarray platforms and
RNA-seq. Sampling time was recovered from GEO sample titles and characteristics for every sample; no
time point was inferred. Two *S. elongatus* PCC 7942 constant-light series (GSE14225, Ito et al.;
GSE18902, Vijayan et al.) from different laboratories on different platforms form the within-strain
benchmark. The *kaiABC*-null arm of GSE14225 is the negative control. *Halothece* sp. PCC 7418
(GSE205640) is reported but excluded from cross-lineage analysis (§3.7).

### 2.2 Rhythm detection

Per gene, log₂ expression was fitted with

  y = Σₛ bₛ·I(series = s) + A·cos(2πt/T) + B·sin(2πt/T),  T = 24 h,

with per-replicate and per-platform intercepts. Significance is an F-test of the harmonic terms
against the series-only model, Benjamini–Hochberg corrected across genes. Phase is the fitted peak.
A 20–28 h period scan was run as a sensitivity check and did not change conclusions.

### 2.3 Orthology

Two independent assignments. (i) KEGG Orthology groups single-copy in the genome annotation and
unique within the dataset. (ii) Reciprocal best hits by Smith–Waterman alignment (pyswrd, BLOSUM62,
E < 10⁻⁵) between RefSeq proteomes. Datasets are quantified against original locus tags while
current proteomes carry modern tags; the RefSeq GFF `old_locus_tag` attribute bridges them exactly.
Neither assignment is tree-based and neither is synteny-validated (§6).

### 2.4 Statistic

Mean |phase difference| between orthologs, minimised over a global shift of −12 to +12 h in 1 h
steps, since studies set their time origin differently. We report both the absolute reduction in
hours and the percentage reduction relative to the null mean; §3.5 shows why the percentage form
alone is misleading.

### 2.5 Nulls

Four, kept distinct:

| null | construction | what it holds fixed |
|---|---|---|
| **uniform-phase** | phases drawn uniformly at random | nothing about the data's phase distribution |
| **time-point permutation** ("Ferrari-like") | permute sample time labels, re-fit rhythms, recompute | the expression values, not the phase distribution |
| **ortholog-identity permutation** (primary) | keep fitted phases, shuffle which ortholog carries which | each dataset's phase distribution exactly |
| **structure-preserving** | circular rotation along the genome; operon-block permutation | the above, plus local phase autocorrelation |

All are minimised over global shifts identically to the observed statistic.

**Operon blocks.** Orthologs are ordered by genomic coordinate in species A. A block break is placed
wherever consecutive orthologs lie on a different contig, a different strand, or more than 2 kb
apart (gene-start to gene-start, from the RefSeq GFF). Blocks are permuted whole, genes within a
block keeping their order and adjacency. No independent operon-prediction database was used, and
the permutation does not constrain where blocks land relative to genomic regions — both are
limitations (§6). For the single-copy ortholog sets analysed, mean block size is ≈1.3 genes, so this
null is close to the gene-level one by construction; its value is as a bound, not as a strong test.

### 2.6 Functional stratification

As a sensitivity analysis, ortholog identity is permuted only within a broad functional class taken
from the KEGG BRITE `ko00001` hierarchy, which is built without reference to any expression
measurement. Level A gives 6–7 classes, level B 22–28. BRITE is not a partition — 38.6%
of KOs appear under more than one level-B class — — 33.3% under more than one level-A class, up to
26 level-B classes for a single KO — so one class per ortholog must be chosen by rule.
First occurrence in the hierarchy file was fixed as the primary rule before any stratified p value
was examined; the whole analysis is then repeated under the last occurrence, an equally arbitrary
rule, as a sensitivity check (§3.8). Neither rule was selected by its outcome, and the sensitivity
run is reported whichever way it came out. Orthologs whose KO carries no BRITE class are kept in an explicit `unclassified`
stratum rather than discarded; on this panel that stratum is empty, because every ortholog entering
the comparison is a KEGG KO that BRITE classifies. This replaces one exchangeability assumption
with a weaker one — orthologs exchangeable *within* a functional class — and that assumption is
still an assumption: broad classes can hide systematic differences in expression level and rhythm
amplitude. Class counts, median class sizes and the number of orthologs in classes smaller than ten
are reported for every pair, together with a sensitivity run that drops those small classes
entirely. The seven pairs examined were fixed before the analysis was run, as the within-strain
benchmark plus every cross-lineage pair involving at least 150 shared orthologs. Every stratified
p value is an empirical tail probability over 1000 permutations, so its Monte Carlo resolution is
0.001 and a reported p of 0.000 means "below 0.001", not zero.

The question this null answers is conditional and differs from the primary one: *is phase
correspondence greater than expected given both the marginal phase distributions and broad
functional composition?* A weaker result under it does not invalidate the primary result; it
localises part of the correspondence to functional composition rather than ortholog identity.

### 2.7 Simulations

Three, all with truth fixed by construction. (i) **Null comparison** (§3.6): sinusoidal expression is
generated for each gene at the sampled time points with noise, phases drawn independently for the
two species from their own concentrated marginals — zero ortholog-specific conservation — and rhythm
detection is run exactly as on real data, so every rejection is a false positive. (ii) **Power**
(§3.4): a fraction *c* of orthologs are assigned conserved timing (phase of B = phase of A + offset
+ Gaussian jitter of SD σ), the rest drawn independently; *c*, σ and n are swept. (iii)
**Design-matched calibration** (§3.6): the same zero-conservation construction run under the
sampling design of each analysed dataset, and separately under each sampling design reported by
Ferrari et al., with inner permutations set to 1000 to match the published procedure. The second set
is **design-matched only**: their empirical phase distributions and expression characteristics are
not used, because that requires their data. It therefore characterises the behaviour of a procedure
under conditions resembling their experiments and says nothing about their result.

Monte Carlo accounting. Outer trials (independent simulated datasets) and inner permutations (the
null distribution within one trial) are reported separately throughout, because they bound different
things: the outer count sets the precision of the reported rate, the inner count sets the resolution
of the p value within each trial. Every rate carries a Wilson 95% interval on the outer count. Each
grid cell draws from a generator seeded from the cell's own coordinates rather than from one stream
advanced across cells, so a cell's result does not depend on how many cells preceded it and an
interrupted run resumes to the same numbers.

---

## 3. Results

### 3.1 The rhythm caller is calibrated

The *kaiABC*-null arm yields **0 of 2490 genes** at q < 0.05. In wild type, *kaiB*/*kaiC* show 3.07/3.07 log₂
rhythms peaking at 11.3/11.4 h (q = 1.3e-09/8.8e-09);
*Trichodesmium kaiC* peaks at 13.8 h. The two amplitudes are equal only after rounding
(*kaiB* Synpcc7942_1217, *kaiC* Synpcc7942_1216; 3.0656 and 3.0673 log₂ unrounded) — they are
co-transcribed from the *kaiBC* operon, so near-identical amplitude is expected rather than a repeated
value. *kaiA* is absent from *Prochlorococcus*, as expected.
Rhythmic fractions are 52–68% in four lineages and 98% in *Prochlorococcus* under light/dark.

### 3.2 Within one strain, phase replicates and set membership does not

Across 2490 shared genes between the two PCC 7942 studies: 1400 and 1324 rhythmic, 966 in both
against 744 expected under independence (Jaccard 0.549). Median |phase difference| is **1.67 h**,
93.8% agreeing within 4 h. This reproduces in a bacterium, in constant light, the dissociation
Brooks et al. reported for mouse liver.

### 3.3 Cross-lineage transfer is an order of magnitude weaker (Figures 1, 2)

| ortholog selection | within strain | between lineages, median | max | pairs significant |
|---|---|---|---|---|
| all one-to-one orthologs, no rhythm threshold | 50.4% | 4.0% | 7.5% | 6 of 10 |
| either dataset q < 0.05 | 61.6% | 4.2% | 8.4% | 6 of 10 |
| both datasets q < 0.05 | 75.0% | 5.0% | 12.3% | 7 of 10 |
| top 200 by combined rank | 82.6% | 6.4% | 13.7% | 6 of 10 |

In absolute terms the within-strain reduction is **4.15 h**; cross-lineage reductions are 0.07–0.55 h
(median by group, §3.5). Reciprocal best hits (1175–1804 orthologs per pair) reproduce the KO-based
values almost exactly:

| pair | RBH orthologs | effect, all orthologs | effect, both q<0.05 | same pair by KO |
|---|---|---|---|---|
| F: 12 pts, 2 h, 1 day | 12 | 12 | 0.010 (0.003-0.029) | 0.177 (0.138-0.224) | 0.433 (0.378-0.490) | 0.547 (0.490-0.602) |
| F: 24 pts, 2 h, 2 days | 24 | 24 | 0.020 (0.009-0.043) | 0.207 (0.165-0.256) | 0.503 (0.447-0.560) | 0.907 (0.868-0.935) |
| F: 7 pts, 4 h (Arab/rice) | 7 | 7 | 0.007 (0.002-0.024) | 0.107 (0.077-0.147) | 0.160 (0.123-0.206) | 0.087 (0.060-0.124) |

**Multiple testing.** Benjamini–Hochberg correction across the ten cross-lineage pairs within each
selection rule, applied as an exploratory adjustment that does not remove the dependence between
pairs sharing a genome, leaves 6–7 of 10 significant at q < 0.05 — the correction does not change the count,
because the significant p values are small and the non-significant ones large. This is reported in
place of the raw count, but the pairs are not independent and §3.5 should be read before treating
"6–7 of 10" as evidence of a group-wide phenomenon.

### 3.4 Power: what the design can and cannot establish

Power at α = 0.05, with 60% of genes peaking in one 6 h window. 63 cells,
400 outer trials and 400 inner permutations each; Wilson 95% intervals in
parentheses:

| true conserved fraction | jitter SD | n=150 | n=500 | n=1500 |
|---|---|---|---|---|
| 0% | 1 h | 0.05 (0.035-0.079) | 0.03 (0.017-0.052) | 0.06 (0.037-0.082) |
| 0% | 2 h | 0.04 (0.021-0.058) | 0.05 (0.035-0.079) | 0.04 (0.023-0.061) |
| 0% | 4 h | 0.07 (0.045-0.094) | 0.06 (0.039-0.085) | 0.03 (0.017-0.052) |
| 10% | 1 h | 0.15 (0.123-0.194) | 0.34 (0.298-0.390) | 0.45 (0.397-0.494) |
| 10% | 2 h | 0.12 (0.092-0.156) | 0.25 (0.210-0.295) | 0.43 (0.387-0.484) |
| 10% | 4 h | 0.08 (0.057-0.111) | 0.15 (0.123-0.194) | 0.32 (0.276-0.367) |
| 20% | 1 h | 0.43 (0.382-0.479) | 0.53 (0.479-0.576) | 0.56 (0.516-0.613) |
| 20% | 2 h | 0.39 (0.341-0.436) | 0.49 (0.444-0.541) | 0.60 (0.551-0.647) |
| 20% | 4 h | 0.25 (0.212-0.297) | 0.39 (0.341-0.436) | 0.45 (0.402-0.499) |
| 40% | 1 h | 0.97 (0.945-0.981) | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) |
| 40% | 2 h | 0.87 (0.831-0.897) | 0.97 (0.951-0.985) | 1.00 (0.990-1.000) |
| 40% | 4 h | 0.53 (0.481-0.578) | 0.72 (0.674-0.762) | 0.74 (0.697-0.783) |
| 60% | 1 h | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) |
| 60% | 2 h | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) |
| 60% | 4 h | 0.98 (0.964-0.991) | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) |
| 80% | 1 h | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) |
| 80% | 2 h | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) |
| 80% | 4 h | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) |
| 100% | 1 h | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) |
| 100% | 2 h | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) |
| 100% | 4 h | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) | 1.00 (0.990-1.000) |

False-positive rate at c = 0 is 0.03–0.07 across every
zero-conservation cell, so the test holds its size. The
design excludes conservation of ≳40% of orthologs when jitter is tight (power 0.87–1.00 at 1–2 h
across all three sample sizes), but that bound weakens substantially under loose jitter
(0.53–0.74 at 4 h). It cannot distinguish 10% from 20% at any setting (power 0.08–0.60). **The mapping from
observed effect to true conserved fraction is a property of this simulation, not a measurement**:
the same ~5% observed effect is compatible with 10–20% of orthologs conserved at 1–2 h jitter and
with a larger fraction conserved more loosely. We therefore state an upper bound and decline to
report a point estimate of the residual.

### 3.5 The residual is unevenly distributed and confounded with phase concentration

The cross-lineage signal is not spread evenly across pairs (both q<0.05, effect %):

| pair | orthologs | effect % | mean phase concentration |
|---|---|---|---|
| *S. elongatus* × *Trichodesmium* | 376 | +12.3 | 0.41 |
| *Synechocystis* × *Nostoc* | 171 | +11.9 | 0.58 |
| *S. elongatus* × *Synechocystis* | 146 | +10.7 | 0.52 |
| *S. elongatus* × *Nostoc* | 338 | +9.5 | 0.51 |
| *Trichodesmium* × *Nostoc* | 375 | +7.0 | 0.47 |
| *S. elongatus* × *Prochlorococcus* | 475 | +2.9 | 0.68 |
| *Prochlorococcus* × *Nostoc* | 447 | +2.9 | 0.74 |
| *Trichodesmium* × *Prochlorococcus* | 488 | -0.5 | 0.64 |
| *Synechocystis* × *Trichodesmium* | 200 | -1.2 | 0.48 |
| *Synechocystis* × *Prochlorococcus* | 219 | -6.2 | 0.75 |

Leave-one-lineage-out, read as a sensitivity analysis over five lineages rather than a test of
which lineage causes the effect: dropping *S. elongatus* lowers it from
5.0% to 1.2%; dropping *Nostoc* to 1.2%; dropping *Prochlorococcus* raises it to 10.1%. In absolute
hours, pairs involving *S. elongatus* gain a median 0.55 h against 0.07 h for pairs without it, so
the dependence is not purely an artifact of the percentage scale.

**But the percentage scale is confounded.** As an exploratory observation on ten non-independent
pairs — not a demonstration of the mechanism — effect size correlates negatively with how concentrated
the two lineages' phase distributions are (Spearman ρ = -0.62 to -0.55 across rules)
and positively with the null mean |phase difference| (ρ = +0.48 to +0.79). Where pile-up is extreme
the null mean is already small and there is little room for a percentage reduction. *Prochlorococcus*
has the most extreme concentration of any dataset (91% of rhythmic transcripts within one 6 h
window) and appears in three of the four weakest pairs.

We therefore cannot separate two explanations: that *Prochlorococcus* — a reduced-genome, *kaiA*-less
lineage — genuinely shares less timing, or that the measure loses sensitivity exactly there.
Excluding it, the median cross-lineage effect is
5.7–10.1% across the four selection rules, with 5 of
6 pairs significant under the both-rhythmic rule. Both
readings are stated; the data do not choose between them.

### 3.6 Three nulls, separated

On simulated expression with conservation fixed at zero: 8 cells, 1000 outer
trials and 200 inner permutations each, 24 time points over 2 days, cosinor detection as on
real data; rejection rate at nominal α = 0.05, from `results/n13_three_nulls_1000.csv`:

| n orthologs | phase concentration | uniform-phase | time-point permutation (95% CI) | ortholog-identity |
|---|---|---|---|---|
| 300 | none (uniform truth) | 0.04 | **0.02** (0.011-0.028) | 0.05 |
| 300 | 40% in a 6 h window | 0.80 | **0.21** (0.184-0.234) | 0.05 |
| 300 | 60% in a 6 h window | 1.00 | **0.55** (0.522-0.584) | 0.05 |
| 300 | 90% in a 6 h window | 1.00 | **0.83** (0.811-0.857) | 0.06 |
| 1000 | none (uniform truth) | 0.05 | **0.01** (0.008-0.022) | 0.04 |
| 1000 | 40% in a 6 h window | 1.00 | **0.48** (0.445-0.507) | 0.05 |
| 1000 | 60% in a 6 h window | 1.00 | **0.76** (0.736-0.788) | 0.05 |
| 1000 | 90% in a 6 h window | 1.00 | **0.92** (0.896-0.931) | 0.04 |

The mechanism is visible in the null means. As phase concentration rises from none to
90%, the uniform null does not move (5.86 →
5.86 h) because it never looks at the data; the time-point-permutation null
adapts partially (5.86 → 4.37 h); the
ortholog-identity null adapts fully (5.86 → 2.85 h).
A null that does not track the marginal cannot distinguish shared pile-up from shared
ortholog-specific timing, which is the whole of the effect being measured.

**Calibration under the designs actually used.** Repeating the zero-conservation construction under
the sampling design of each analysed dataset. Across both panels below: 22 cells,
300 outer trials and 1000 inner permutations per cell, 300 orthologs,
rejection rate at nominal α = 0.05 with Wilson 95% intervals in parentheses. Phases are drawn
independently for the two datasets, so the true conserved fraction is zero by construction and
every rejection is a false positive; concentration is the fraction of genes whose true phase falls
in one 6 h window. The null hypothesis tested in each cell is that the two phase vectors are
unrelated given the design, and the test is the one-sided empirical tail probability of the
observed mean |phase difference| under the stated null:

| sampling design | time points | samples | uniform truth | 60% in 6 h |
|---|---|---|---|---|
| ours: MED4 | 12 | 36 | 0.033 (0.018-0.060) | 0.033 (0.018-0.060) |
| ours: NPUN | 10 | 30 | 0.060 (0.038-0.093) | 0.033 (0.018-0.060) |
| ours: SE7942_ito_WT | 13 | 24 | 0.043 (0.025-0.073) | 0.040 (0.023-0.069) |
| ours: SYN6803 | 24 | 24 | 0.047 (0.028-0.077) | 0.063 (0.041-0.097) |
| ours: TERY | 9 | 27 | 0.040 (0.023-0.069) | 0.047 (0.028-0.077) |

The primary null holds its nominal size under every design this manuscript relies on.

Under the three sampling designs reported by Ferrari et al., with the time-point-permutation null at
the published 1000 permutations:

| sampling design | time points | samples | uniform truth | 40% in 6 h | 60% in 6 h | 90% in 6 h |
|---|---|---|---|---|---|---|
| F: 12 pts, 2 h, 1 day | 12 | 12 | 0.010 (0.003-0.029) | 0.177 (0.138-0.224) | 0.433 (0.378-0.490) | 0.547 (0.490-0.602) |
| F: 24 pts, 2 h, 2 days | 24 | 24 | 0.020 (0.009-0.043) | 0.207 (0.165-0.256) | 0.503 (0.447-0.560) | 0.907 (0.868-0.935) |
| F: 7 pts, 4 h (Arab/rice) | 7 | 7 | 0.007 (0.002-0.024) | 0.107 (0.077-0.147) | 0.160 (0.123-0.206) | 0.087 (0.060-0.124) |

**This is design-matched, not data-matched.** Phase concentration is swept rather than taken from
their data, so these rows describe how the procedure behaves under conditions resembling their
experiments. They are a property of a sampling design and a null, not of anybody's result. No rate
reported here is attributable to the published analysis, which used different organisms, different
data and a different rhythm caller; attaching any of these numbers to it would be unsupported.

**The coarse design is not monotonic in phase concentration, and must not be reduced to one
inflation factor.** The 7-point, 4 h design (row `F: 7 pts, 4 h (Arab/rice)` above) gives a
false-positive rate of
0.007 (0.002-0.024) at uniform phases, 0.107 (0.077-0.147) at 40%
concentration, 0.160 (0.123-0.206) at 60%, and 0.087 (0.060-0.124) at 90%.
The rate therefore peaks at intermediate concentration and falls again at high concentration, and
the 60% and 90% intervals do not overlap, so this is not Monte Carlo noise. Each design and each
concentration must be read on its own row; no single correction factor describes the family, and
no null model is generally superior across it.

The mechanism we would propose is that two quantities move together, and the statistic is their
difference. Seven points 4 h apart resolve phase onto a coarse grid. As the phase distribution
concentrates, the observed mean |phase difference| falls, but so does the mean of the
time-point-permutation null, because refits after permuting time labels also land on the same few
grid positions. The reported rate reflects the gap between them, which can widen and then narrow.
That the null's own location moves with concentration is directly measured, though on the 24-point
design rather than this one: the time-point null's mean |phase difference| falls from
5.86 h at uniform phases to 4.37 h at 90%
concentration, in the three-null table above. **We have not measured the observed and null locations separately within the
7-point design**, so the mechanism is an explanation consistent with the numbers rather
than one this simulation establishes; separating them would require recording both locations per
cell, which these runs do not.

**This also corrects v0.1, which described the published comparison's null as expecting uniform phases.**
It does not. A time-point-permutation null is substantially protected relative to a uniform null,
and it is the appropriate object of criticism. It is nonetheless inflated — a false-positive rate of
0.21–0.92 against a nominal 0.05, and like the uniform null it worsens with more
orthologs (0.21→0.48 at 40% concentration,
0.83→0.92 at 90%, from n=300 to n=1000) — across the entire range of phase concentration observed in real
data (0.32–0.91 here).

### 3.7 Structure, marginals, and the limits of the panel

Genomically adjacent same-strand genes within 2 kb share phase (ratio of adjacent to random mean |Δphase| 0.41–0.73; 27–59% of adjacent pairs within 1 h), so genes are not exchangeable
genome-wide. Repeating the tests under circular rotation and operon-block permutation changes
effects and p values negligibly (11.2%/p=0.003 → 11.1%/p=0.008 → 11.1%/p=0.007), because the
single-copy ortholog subset is sparse along the genome (§2.5).

Each lineage concentrates transcripts at its own time of day, from 32% to 91% in the busiest 6 h
window, with non-coincident windows (Figure 3).

Only two lineages have genome-wide free-running data. The *Halothece* series gives 6 time points
over 30 h per condition, 0.8–4.0% rhythmic genes, and a best-fit period
pinned at the 28 h scan boundary — underpowered. Effectively one cyanobacterial lineage has adequate constant-light data.

### 3.8 Conditioning on functional composition

Permuting ortholog identity only within broad functional classes asks whether correspondence is
carried by the orthologs themselves or by the functional composition of the sets being compared.

At BRITE level A (6–7 classes):

| pair | orthologs | classes | median class size | free null | stratified | n in classes ≥10 | stratified, classes ≥10 |
|---|---|---|---|---|---|---|---|
| **within strain, 7942 ×2** | 480 | 7 | 57 | 75.0% | 73.5% (p=0.000) | 475 | 73.5% (p=0.000) |
| *S. elongatus* × *Synechocystis* | 146 | 6 | 18 | 10.7% | 9.7% (p=0.010) | 137 | 11.6% (p=0.003) |
| *S. elongatus* × *Trichodesmium* | 376 | 7 | 40 | 12.3% | 12.2% (p=0.000) | 365 | 12.0% (p=0.000) |
| *S. elongatus* × *Prochlorococcus* | 475 | 7 | 42 | 2.9% | 1.4% (p=0.129) | 469 | 1.3% (p=0.155) |
| *S. elongatus* × *Nostoc* | 338 | 7 | 40 | 9.5% | 7.2% (p=0.008) | 330 | 6.1% (p=0.030) |
| *Trichodesmium* × *Nostoc* | 375 | 7 | 41 | 7.0% | 6.1% (p=0.007) | 371 | 5.9% (p=0.007) |
| *Synechocystis* × *Nostoc* | 171 | 6 | 22 | 11.9% | 10.6% (p=0.005) | 163 | 10.9% (p=0.012) |

At BRITE level B (22–28 classes):

| pair | orthologs | classes | median class size | free null | stratified | n in classes ≥10 | stratified, classes ≥10 |
|---|---|---|---|---|---|---|---|
| **within strain, 7942 ×2** | 480 | 28 | 12 | 75.0% | 71.9% (p=0.000) | 432 | 72.1% (p=0.000) |
| *S. elongatus* × *Synechocystis* | 146 | 23 | 4 | 10.7% | 8.1% (p=0.017) | 75 | 8.0% (p=0.068) |
| *S. elongatus* × *Trichodesmium* | 376 | 27 | 9 | 12.3% | 11.0% (p=0.000) | 319 | 10.2% (p=0.000) |
| *S. elongatus* × *Prochlorococcus* | 475 | 26 | 10 | 2.9% | 0.8% (p=0.246) | 422 | 1.1% (p=0.202) |
| *S. elongatus* × *Nostoc* | 338 | 23 | 9 | 9.5% | 6.2% (p=0.028) | 282 | 2.9% (p=0.212) |
| *Trichodesmium* × *Nostoc* | 375 | 24 | 11 | 7.0% | 6.6% (p=0.002) | 319 | 4.7% (p=0.031) |
| *Synechocystis* × *Nostoc* | 171 | 22 | 6 | 11.9% | 9.6% (p=0.014) | 93 | 5.2% (p=0.195) |

Three things follow. First, the within-strain benchmark is essentially unaffected at every level
(73.5%, 71.9%, 72.1% with small classes dropped), so the
reference scale is not an artifact of functional composition. Second, at coarse classes
5 of 6 cross-lineage pairs remain significant, and at fine classes
5 of 6 do; with small strata removed only 2 of
6 do. Third, which pairs those are is not stable, so the aggregate count is the only
level at which we read this analysis; the next paragraph gives the reason and the evidence.

The whole stratified comparison is **exploratory**. It conditions on an annotation hierarchy built
for a different purpose, the pairs are not independent, and the quantity it reports is a count over
seven dependent comparisons, not a test with a controlled error rate.

**Per-pair stratified p values depend on an arbitrary annotation choice.** BRITE is not a partition:
33.3% of KOs sit under more than one level-A category and 38.6% under more
than one level-B category, up to 26 classes for a single KO. Stratification needs one
class per ortholog, so a rule is required, and the reported analysis takes each KO's first
occurrence in the hierarchy file. Repeating everything with the last occurrence instead — an
equally defensible and equally arbitrary rule — leaves the count of significant cross-lineage pairs
almost unchanged (5 → 6 at level A, 5 →
5 at level B) but moves 1 pair across α = 0.05 at level A and
2 at level B. The pair that moves at both levels is *S. elongatus* ×
*Prochlorococcus*: at level A, p = 0.129 under the first-occurrence rule and
0.009 under the last, with the effect moving from 1.4% to
3.2%; at level B, 0.246 and 0.000, effect
0.8% to 3.7%. The pairs that move are, at level A,
*S. elongatus* × *Prochlorococcus* (p 0.129 vs 0.009); at level B, *S. elongatus* × *Prochlorococcus* (p 0.246 vs 0.000); *S. elongatus* × *Nostoc* (p 0.028 vs 0.214). Its effect is small enough that
bookkeeping over multiply-annotated genes decides its significance. We therefore report the
aggregate behaviour of the stratified null and **do not** interpret any individual pair's stratified
p value, including the *Prochlorococcus* one. This supersedes a claim in an earlier draft that the
pair failed under every stratified variant.

What "stable under the assignment rule" means here is narrow and worth stating exactly, because the
word can be stretched to cover three different things. It means: the **count** of cross-lineage
pairs significant at α = 0.05 changes by at most one when the rule is swapped, and the **sign** of
every pair's effect is unchanged. It does **not** mean the identity of the significant pairs is
stable — 1 and 2 pairs cross α — nor that effect **sizes** are
stable, which for small-effect pairs they are not. The first-occurrence rule was fixed as primary
before any stratified p value was examined, and was not chosen for the result it gives; under the
alternative rule the aggregate count would support the same reading.

**The attenuation at fine strata is confounded with power.** Dropping classes smaller than ten
removes 53–78 orthologs per pair and in two cases roughly halves the
comparison (*S. elongatus* × *Synechocystis* 146 → 75;
*Synechocystis* × *Nostoc* 171 → 93), so loss of significance cannot
be attributed to functional composition alone.

Three readings remain open, and this panel does not separate them. (i) The correspondence reflects
functional composition, and conditioning on it removes a confound. (ii) The correspondence is
ortholog-specific but weak, surviving coarse conditioning and falling below detection under fine
conditioning. (iii) The test loses power under fine conditioning, because small strata admit few
distinguishable permutations and the permutable set shrinks, so significance is lost for a reason
unrelated to functional composition. **Loss of significance after stratification is therefore
evidence of sensitivity to functional conditioning, not evidence that ortholog-specific similarity
is absent.** The defensible statement is that the ortholog-identity-specific component is robust to
coarse functional conditioning, and is not established independently of functional composition at
fine conditioning on this panel.

### 3.9 A candidate phase-conserved set, reported without an inferential test

Aligning lineages by best global shift and ranking orthologs by phase dispersion yields 33
tightly co-phased orthologs — dispersion below 2.74 h, rhythmic in at least
4 of the five lineages — with biosynthetic annotation (*fabI*, *fabZ*, *dxr*,
*ispH*, *metK*, *ahcY*, *sat*, *glgC*, *pdhB*, *yidC*), superficially a pre-dawn anabolic module.
**The identified genes constitute an exploratory candidate set; their enrichment has not been
established by a valid selection-adjusted null test.** We attach no significance claim to it.

The reason is that we have no valid test of it. The natural control shuffles ortholog identity
independently within each lineage and recounts the orthologs falling below the same dispersion
cutoff. But both the cutoff and the set are reached through a global alignment shift fitted to the
observed ortholog pairings, and that shift is held fixed across permutations. The null therefore
does not reproduce the selection procedure that produced the observed set, and a null that omits
the selection step cannot calibrate a statistic that depends on it: it is anti-conservative by an
amount this design cannot quantify. **We therefore report no p value and no fold enrichment for
this set.** The counts from the fixed-alignment shuffle are recorded in
`results/n5_core_sd_null.csv` at both membership criteria, as a diagnostic and not as a test.

Establishing enrichment would require refitting the alignment inside every permutation, or fixing
the shift by a rule independent of these data. Neither is done here (§6, gap 13). Pathway
enrichment of the candidate set reaches q < 0.05 only for a generic superset term, so its
annotation does not single out a pathway either.

**This corrects v0.2**, which reported a control computed at a different membership criterion
(rhythmic in at least three lineages) from the core it described, and drew the conclusion "no
identifiable conserved gene set" from that mismatched comparison. The mismatch is fixed. The
conclusion has not been replaced by its opposite, because the corrected comparison rests on the
same invalid null.

---

## 4. Discussion

**What the data support.** Diel transcriptional programs transfer between cyanobacterial lineages at
a small fraction — on the order of a tenth — of what transfers between two independent studies of
one strain. Under tight phase-jitter assumptions (1–2 h) the design has high power to detect timing
conservation involving 40% of orthologs, and none is seen; that exclusion weakens substantially at
4 h jitter. The residual is small,
unevenly distributed across pairs, and partly confounded with marginal phase concentration.

**What the data do not support.** We cannot attribute the difference to rewiring of clock output.
Six of eight time courses are light/dark, where a gene may peak in the morning because it responds
to light rather than because the oscillator gates it; only *S. elongatus* has adequate free-running
data. Rewiring of clock-to-network connections is one hypothesis; divergent direct light responses,
divergent metabolic constraints, and differences in growth condition between the source studies are
others, and this design does not discriminate among them. The single-gene observation of Kennedy et
al. [2025], that a homologue of the canonical circadian input kinase is arrhythmic in *Nostoc
punctiforme*, is consistent with rewiring but does not establish it either.

**This paper answers a question about diel transcriptional phase, and cannot settle the broader
question of whether endogenous circadian regulation is conserved.** Diel expression under light and
dark confounds oscillator-gated timing with direct light response, and only one lineage here has
adequate constant-light data, so nothing reported above either establishes or excludes conserved
downstream circadian regulation. Settling that needs free-running time courses in more than one
lineage; it is a limit of what these data can address, not a shortcoming of the analysis applied
to them.

Two narrower cautions follow the same logic. The stability of the aggregate stratified result under
the BRITE assignment rule is a statement about a count and a sign, not biological robustness: the
identity of the significant pairs and the sizes of their effects are **not** stable, which is why
no individual pair's stratified p value is interpreted. And the phase-conserved orthologs of
§3.9 are an exploratory candidate set with no valid test behind them, offered for targeted
experiments rather than as a finding.

Generalisation to Cyanobacteria as a phylum is limited by five genomes, by their non-independence,
and by the absence of free-running data outside one strain.

**The methodological result.** Any comparison of event timing between groups with concentrated
marginal distributions faces the hazard quantified in §3.6, and the correction is inexpensive:
preserve the marginals, permute only the pairing. Our simulation establishes that a
time-point-permutation null is inflated under realistic phase concentration. It does **not**
establish that Ferrari et al.'s conclusion is wrong: their statistic, orthologs, organisms and
phase-concentration levels all differ from the simulation, and when both nulls were run on our data
they gave similar verdicts. The claim is that re-analysis under a marginal-preserving null is
warranted, and that a re-analysis finding the original conclusion intact would be as informative as
one overturning it.

**Priority.** A documented search (8 PubMed queries, 265 unique records
screened, 31 adjudicated; all 34 records citing Ferrari et al. screened at
title and abstract, of which 0 mention a cyanobacterial taxon and
0 mention a null model) found no study meeting the primary criterion. The strategy
and its five limitations are in `LITERATURE-SEARCH.md`. We state only that to our knowledge this has
not been tested within cyanobacterial lineages, and rest the contribution on the comparison
framework rather than on priority.

**Priority next steps.** (i) Re-run the Archaeplastida comparison on the original data under a
marginal-preserving null (ArrayExpress E-MTAB-7188; Supplementary Data 2–11 of Ferrari et al.);
reproduce the published analysis first and change one component at a time. (ii) Generate or locate
free-running time courses outside PCC 7942 — §3.7 identifies this as the binding constraint on the
biological question. (iii) Extend the panel to resolve 10% from 20%.

---

## 5. Figures, tables and supplementary material

**Figure 1** Effect size per lineage pair under four selection rules, within-strain benchmark
highlighted. `figures/fig1_effect_size.png`
**Figure 2** Cumulative distribution of |phase difference| after best global shift. `figures/fig2_phase_difference_ecdf.png`
**Figure 3** Marginal phase distributions per dataset. `figures/fig3_phase_marginals.png`
**Table 1** Dataset panel. `results/dataset_summary.csv`

Supplementary tables in `results/`: per-gene rhythm calls; clock-gene control; within-species
reproducibility; KO conservation; rank and amplitude sweeps; pairwise phase tests with and without
offset optimisation; two-null and three-null comparisons; threshold sensitivity; RBH comparison;
null simulation; power analysis and power grid; independence check; block nulls; BH-corrected pairs
and leave-one-out; effect-versus-pile-up; conserved-core analysis and pathway enrichment.

Reproduction: `make check && make all` with `requirements.txt` pinned; inputs download to
`$DIEL_DATA_DIR` (default `./data`).

---

## 6. Known gaps

1. The attenuation of cross-lineage correspondence under fine functional strata (§3.8) is
   confounded with the loss of orthologs that stratification causes; the two cannot be separated on
   this panel.
2. The functional-stratified null assumes exchangeability *within* a BRITE class, which broad
   classes may violate through systematic differences in expression level or rhythm amplitude.
3. BRITE assigns 38.6% of KOs to several classes, and which class a stratification uses
   is an arbitrary choice that decides the significance of the weakest pairs (§3.8). Only the
   aggregate behaviour of the stratified null is reported for that reason; a pair-level stratified
   claim would need an annotation source that is a genuine partition.
4. Design-matched simulation of the published sampling schemes uses swept phase concentrations, not
   their empirical distributions; no statement about Ferrari et al.'s actual result follows from it.
5. Conservation restricted to a single functional module rather than spread across orthologs is not
   simulated and therefore not excluded by the power surface.

6. No phylogenetic correction; ten pairs from five genomes, and §3.5 shows the result is sensitive
   to which lineages are included. Five taxa are too few for phylogenetic regression, so the pairs
   are reported as interdependent descriptive comparisons.
7. Orthology is RBH and single-copy KO — neither tree-based nor synteny-validated.
8. Operon blocks are defined by distance and strand only, without an independent operon-prediction
   database, and block permutation does not constrain genomic placement.
9. Power simulation assumes Gaussian phase jitter and uniform conservation across functional
   categories; conservation restricted to particular pathways was not simulated.
10. Phase-estimation uncertainty is not propagated into the cross-lineage statistic.
11. Figures are analysis-grade; no reference list in journal format; no author list, funding or
   ethics statements.
12. Prior-art search covered PubMed, one extended web search, and the 34 papers citing Ferrari et al.
   2019 (16 screened by title for timing/ortholog terms; none compares cyanobacterial lineages or
   uses a marginal-preserving null). bioRxiv full text was not searched.
13. There is no valid test of the candidate phase-conserved set (§3.9). The set is selected
   through a global alignment shift fitted to the observed ortholog pairings, and the available
   control holds that shift fixed across permutations, so it does not reproduce the selection
   procedure and is anti-conservative by an amount the design does not quantify. Testing the set
   would require refitting the alignment inside every permutation, or fixing the shift by a rule
   independent of these data. Neither is done here, and no p value or enrichment is reported for
   the set in either direction.
14. The automated check that forbids hand-typed statistics in the Abstract and Results treats any
   decimal, and any integer of three or more digits, as a statistic. One- and two-digit integers
   are still read as prose and are not caught, so small counts in those sections rest on the
   semantic assertions and on review rather than on that check.

## 7. Changes from v0.2

- **Every statistic is now generated** from `results/` via `src/render_manuscript.py`;
  `src/check_manuscript.py` fails on drift or on any hand-typed statistic in the Abstract or
  Results. Switching this on exposed four errors in v0.2 (series counted as arms; kaiB values
  attributed to kaiC; a *Trichodesmium* phase from the wrong gene row; an RBH range that included
  the within-strain identity case).
- **Semantic assertions added** (`src/validate_semantics.py`, `make validate`): clock-gene labels
  must match their KEGG KO, arm and series counts must be derived from the registry, accessions must
  resolve to matrices that declare them, phases must lie in [0,24), the negative control must stay
  negative, and the manuscript must not contain claims stronger than the results.
- **The fast cosinor rewrite is validated** against the reference implementation
  (`src/test_cosinor_equivalence.py`): amplitude, phase, p and q agree to 0.0e+00.
- **Functional stratification added** (§2.6, §3.8) at two BRITE levels with small-stratum sensitivity.
- **Power surface completed** to 10/20/40/60/80/100% conserved with Monte Carlo intervals.
- **Null calibration under the actual sampling designs** of each dataset, and separately under the
  three published Ferrari designs at 1000 inner permutations, labelled design-matched only.
- **Literature search documented** in `LITERATURE-SEARCH.md`; Lehmann et al. 2013 added as prior
  methodological work on cyanobacterial time-series analysis.
- **§3.9 withdrawn as an inferential result.** v0.2 concluded that no conserved gene set was
  supported, from a control computed at a different membership criterion from the core it
  described. Correcting the mismatch made the comparison significant; that comparison is also
  invalid, because the set is selected through an alignment shift fitted to the observed ortholog
  pairings which the null then holds fixed. No p value or enrichment is reported for the set in
  either direction, and the registry keys for them were removed so they cannot reappear (gap 13).
- **Per-pair stratified p values withdrawn.** BRITE is not a partition (38.6% of KOs
  carry several level-B classes), and the assignment rule moves pairs across α. Only aggregate
  behaviour is interpreted, "stable" is defined explicitly in §3.8, and both the primary and the
  alternative rule are reported with a statement that neither was chosen by its outcome.
- **Stale statistics in §3.4 corrected.** The power prose quoted a 40-trial pilot while the table
  beside it reported the 400-trial surface; the hand-typed literals disagreed with the Abstract's
  templated values for the same quantities. The numeric-literal check that should have caught this
  had exempted every decimal as a "section number" and is repaired.
- **Provenance and coverage are now reported artifacts**: `results/NUMBER_PROVENANCE.csv` maps each
  printed number to its source file, producing target and a test proving the manuscript depends on
  it; `results/ASSERTION_COVERAGE.csv` maps assertions to the mutation tests that prove they fire,
  and reports the two counts separately rather than implying one test per assertion.

## 8. Changes from v0.1

- Primary claim narrowed from circadian clock output to **diel transcriptional programs**; the three
  regulatory concepts are now separated in §1 and the rewiring interpretation demoted to a
  hypothesis in §4.
- Abstract corrected: the published comparison's null is a **time-point permutation**, not a
  uniform-phase null; §3.6 now measures all three separately. "Cross-kingdom" replaced by
  "Archaeplastida".
- Introduction no longer asserts a "universal prokaryotic clock"; the three distinct claims are
  separated.
- Brooks et al. described accurately: technical, biological and statistical contributions, not
  thresholding alone.
- Added: BH correction across pairs; leave-one-lineage-out; the *S. elongatus* dependence; the
  effect-size/phase-concentration confound; absolute as well as percentage effects.
- Power analysis expanded over sample size, jitter and conservation level; the 10–20% figure is
  explicitly labelled a simulation property, not a measurement.

**Correction applied after v0.2 was first written.** Several power values in §3.4 and the upper end
of the time-point-permutation false-positive range had been transcribed from a partially completed
simulation run. Both tables are now generated directly from the result CSVs. The corrected figures
change two statements: power at n=150 for 40% conservation is 0.88–0.93 under 1–2 h jitter rather
than the 0.53–0.60 reported, so the exclusion bound is stronger than stated except under loose
jitter; and the time-point-permutation null's worst false-positive rate is 0.88, not 0.76.
- Operon-block construction documented.
- `requirements.txt` pinned; `Makefile` added as a single entry point.

---

## 9. Declarations

**Author.** Melad Naieem Fathy, independent researcher. ORCID
[0009-0002-4553-8894](https://orcid.org/0009-0002-4553-8894).

**Use of AI tools.** The author found the idea, critiqued it, searched the literature, reviewed
the work, and analysed and interpreted the results. A generative AI assistant also did some of this
work, including writing code. The author takes full responsibility for the content.

**Data availability.** Every input is a public GEO series or NCBI record; no credentials are
required. `src/datasets.py` is the registry of series, arms and lineages, and `src/fetch_data.py`
downloads all of them. The derived result tables are in `results/` in the public release, and
`results/PROVENANCE.json` records a checksum for every input and every output.

**Code availability.** All analysis code is in `src/`, under the MIT licence, with a single entry
point (`make all`). Text, figures and data tables are under CC BY 4.0.

**Reproducibility.** `RELEASE-v0.3.md` lists every released file with a SHA-256.
`CLOSURE-REPORT.md` gives three separate verdicts — content reproducibility, byte-level
reproducibility, and scientific validity — and states what each does and does not establish.

**Pre-registration.** None. This study is exploratory throughout; several central analysis choices
were made with the data in view, and `prereg/NO-PREREGISTRATION.md` lists which ones and when.

**Competing interests.** None.

**Funding.** None. The work was self-funded.
