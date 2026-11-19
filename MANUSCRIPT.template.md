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
lineages has not been measured. We reanalysed {{n_series}} public diel and circadian transcriptome GEO series
({{n_datasets}} analysed arms) from {{n_lineages}} lineages and four platforms, recovering sampling time from deposited sample
metadata in every case. Rhythm detection was validated on a *kaiABC*-null negative control ({{kainull_rhythmic}}
rhythmic genes of {{kainull_genes}} at FDR < 0.05) and on *kaiBC* in wild type ({{kaiB_amp}} log₂ peak-to-trough,
dusk-phased at {{kaiB_phase}}/{{kaiC_phase}} h, q = {{kaiB_q}}/{{kaiC_q}}).

We calibrate cross-lineage comparison against an internal benchmark that cross-species studies have
lacked: two independent studies of the same strain, from different laboratories on different
platforms. Within that strain, ortholog peak times agree to a median of {{within_median_dphase}} h, and the mean
|phase difference| falls {{within_abs_gain_h}} h ({{within_both}}%) below a null preserving each dataset's phase distribution.
Across lineages the same statistic falls {{gain_without_SE_h}}–{{gain_with_SE_h}} h (median {{cross_both_median}}%, maximum {{cross_both_max}}%). The gap holds under
four ortholog-selection rules including none, under two orthology assignments (single-copy KEGG
Orthology; Smith–Waterman reciprocal best hits, {{rbh_min}}–{{rbh_max}} orthologs per pair), and under nulls
preserving operon structure.

Power analysis across sample sizes (150–1500 orthologs), conservation levels and jitter widths shows the
design detects conservation of 40% of orthologs with power {{power40_tight_min}}–{{power40_tight_max}} under tight phase jitter
(1–2 h) but only {{power40_loose_min}}–{{power40_loose_max}} under loose jitter (4 h), and reaches only {{power_small_min}}–{{power_small_max}} at 10–20%. Under tight-jitter assumptions we therefore exclude timing conservation involving 40% of
orthologs; that exclusion is substantially weaker under 4 h jitter, and we cannot resolve the
residual. We do not claim to exclude a majority of the program: 40% is not a majority, and the
bound is conditional on the simulation's assumptions. The residual is also unevenly distributed: it is carried
disproportionately by pairs involving *S. elongatus* (median {{drop_MED4}}% versus {{drop_SE}}% for pairs without it)
and is not detected in pairs involving *Prochlorococcus* — not detected, not shown to be zero —
and effect size correlates negatively, in an exploratory analysis over dependent pairs, with how
concentrated the marginal phase distributions are (ρ = {{pileup_rho}}), so sensitivity and biology
are not separable here.

Conditioning additionally on broad functional composition — permuting ortholog identity only within
KEGG BRITE classes — leaves the within-strain benchmark at {{stratA_within}}% and keeps
{{stratA_sig}} of {{stratA_n}} cross-lineage pairs significant at the coarse class level, but only
{{stratB_sig_big}} of {{stratB_n}} under fine classes with small strata removed, where the
permutable set also shrinks substantially. Individual pairs' stratified p values are not
interpretable: {{multi_pct_B}}% of KEGG Orthology groups carry more than one BRITE class, and
resolving that differently moves {{stratB_n_flip}} of {{stratB_n}} pairs across α = 0.05 without
changing the aggregate. Phase correspondence alone is therefore not sufficient evidence of a
conserved regulatory program.

Separately, we compare three nulls on simulated expression with conservation fixed at zero. A null
expecting uniform phases has a false-positive rate of {{fpr_uniform_min}}–{{fpr_uniform_max}} under realistic phase concentration;
a time-point-permutation null of the kind used in a published Archaeplastida comparison is partially
but not fully protected ({{fpr_timepoint_min}}–{{fpr_timepoint_max}}); an ortholog-identity-permutation null preserving the marginals
holds {{fpr_marginal_min}}–{{fpr_marginal_max}}. Every dataset examined here lies in the regime where the first two are inflated.

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

{{n_series}} GEO series, {{n_datasets}} analysed arms, from {{n_lineages}} lineages (Table 1), 2009–2025, across three microarray platforms and
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
limitations (§6). For the single-copy ortholog sets analysed, mean block size is ≈{{mean_block_size}} genes, so this
null is close to the gene-level one by construction; its value is as a bound, not as a strong test.

### 2.6 Functional stratification

As a sensitivity analysis, ortholog identity is permuted only within a broad functional class taken
from the KEGG BRITE `ko00001` hierarchy, which is built without reference to any expression
measurement. Level A gives 6–7 classes, level B 22–28. BRITE is not a partition — {{multi_pct_B}}%
of KOs appear under more than one level-B class — — {{multi_pct_A}}% under more than one level-A class, up to
{{multi_max_B}} level-B classes for a single KO — so one class per ortholog must be chosen by rule.
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

The *kaiABC*-null arm yields **{{kainull_rhythmic}} of {{kainull_genes}} genes** at q < 0.05. In wild type, *kaiB*/*kaiC* show {{kaiB_amp}}/{{kaiC_amp}} log₂
rhythms peaking at {{kaiB_phase}}/{{kaiC_phase}} h (q = {{kaiB_q}}/{{kaiC_q}});
*Trichodesmium kaiC* peaks at {{kaiC_tery_phase}} h. The two amplitudes are equal only after rounding
(*kaiB* {{kaiB_locus}}, *kaiC* {{kaiC_locus}}; {{kaiB_amp4}} and {{kaiC_amp4}} log₂ unrounded) — they are
co-transcribed from the *kaiBC* operon, so near-identical amplitude is expected rather than a repeated
value. *kaiA* is absent from *Prochlorococcus*, as expected.
Rhythmic fractions are 52–68% in four lineages and 98% in *Prochlorococcus* under light/dark.

### 3.2 Within one strain, phase replicates and set membership does not

Across {{within_genes}} shared genes between the two PCC 7942 studies: {{within_rhy_1}} and {{within_rhy_2}} rhythmic, {{within_rhy_both}} in both
against {{within_rhy_exp}} expected under independence (Jaccard {{within_jaccard}}). Median |phase difference| is **{{within_median_dphase}} h**,
{{within_pct_4h}}% agreeing within 4 h. This reproduces in a bacterium, in constant light, the dissociation
Brooks et al. reported for mouse liver.

### 3.3 Cross-lineage transfer is an order of magnitude weaker (Figures 1, 2)

{{table:effect_table}}

In absolute terms the within-strain reduction is **{{within_abs_gain_h}} h**; cross-lineage reductions are {{gain_without_SE_h}}–{{gain_with_SE_h}} h
(median by group, §3.5). Reciprocal best hits ({{rbh_min}}–{{rbh_max}} orthologs per pair) reproduce the KO-based
values almost exactly:

{{table:rbh_table}}

**Multiple testing.** Benjamini–Hochberg correction across the ten cross-lineage pairs within each
selection rule, applied as an exploratory adjustment that does not remove the dependence between
pairs sharing a genome, leaves {{cross_all_sig}}–{{cross_both_sig}} of 10 significant at q < 0.05 — the correction does not change the count,
because the significant p values are small and the non-significant ones large. This is reported in
place of the raw count, but the pairs are not independent and §3.5 should be read before treating
"{{cross_all_sig}}–{{cross_both_sig}} of 10" as evidence of a group-wide phenomenon.

### 3.4 Power: what the design can and cannot establish

Power at α = 0.05, with 60% of genes peaking in one 6 h window. {{power_cells}} cells,
{{power_outer}} outer trials and {{power_inner}} inner permutations each; Wilson 95% intervals in
parentheses:

{{table:power_surface_table}}

False-positive rate at c = 0 is {{power_size_min}}–{{power_size_max}} across every
zero-conservation cell, so the test holds its size. The
design excludes conservation of ≳40% of orthologs when jitter is tight (power {{power40_tight_min}}–{{power40_tight_max}} at 1–2 h
across all three sample sizes), but that bound weakens substantially under loose jitter
({{power40_loose_min}}–{{power40_loose_max}} at 4 h). It cannot distinguish 10% from 20% at any setting (power {{power_small_min}}–{{power_small_max}}). **The mapping from
observed effect to true conserved fraction is a property of this simulation, not a measurement**:
the same ~5% observed effect is compatible with 10–20% of orthologs conserved at 1–2 h jitter and
with a larger fraction conserved more loosely. We therefore state an upper bound and decline to
report a point estimate of the residual.

### 3.5 The residual is unevenly distributed and confounded with phase concentration

The cross-lineage signal is not spread evenly across pairs (both q<0.05, effect %):

{{table:pair_table}}

Leave-one-lineage-out, read as a sensitivity analysis over five lineages rather than a test of
which lineage causes the effect: dropping *S. elongatus* lowers it from
{{loo_all}}% to {{drop_SE}}%; dropping *Nostoc* to {{drop_NPUN}}%; dropping *Prochlorococcus* raises it to {{drop_MED4}}%. In absolute
hours, pairs involving *S. elongatus* gain a median {{gain_with_SE_h}} h against {{gain_without_SE_h}} h for pairs without it, so
the dependence is not purely an artifact of the percentage scale.

**But the percentage scale is confounded.** As an exploratory observation on ten non-independent
pairs — not a demonstration of the mechanism — effect size correlates negatively with how concentrated
the two lineages' phase distributions are (Spearman ρ = {{pileup_rho}} across rules)
and positively with the null mean |phase difference| (ρ = {{nullmean_rho}}). Where pile-up is extreme
the null mean is already small and there is little room for a percentage reduction. *Prochlorococcus*
has the most extreme concentration of any dataset ({{pileup_max_pct}}% of rhythmic transcripts within one 6 h
window) and appears in three of the four weakest pairs.

We therefore cannot separate two explanations: that *Prochlorococcus* — a reduced-genome, *kaiA*-less
lineage — genuinely shares less timing, or that the measure loses sensitivity exactly there.
Excluding it, the median cross-lineage effect is
{{exmed4_median_min}}–{{exmed4_median_max}}% across the four selection rules, with {{exmed4_nsig}} of
{{exmed4_npairs}} pairs significant under the both-rhythmic rule. Both
readings are stated; the data do not choose between them.

### 3.6 Three nulls, separated

On simulated expression with conservation fixed at zero: {{tn_cells}} cells, {{tn_outer}} outer
trials and {{tn_inner}} inner permutations each, 24 time points over 2 days, cosinor detection as on
real data; rejection rate at nominal α = 0.05, from `results/{{tn_source}}`:

{{table:three_null_table}}

The mechanism is visible in the null means. As phase concentration rises from none to
{{pileup_sim_max}}, the uniform null does not move ({{nullmean_uniform_flat}} →
{{nullmean_uniform_conc}} h) because it never looks at the data; the time-point-permutation null
adapts partially ({{nullmean_timepoint_flat}} → {{nullmean_timepoint_conc}} h); the
ortholog-identity null adapts fully ({{nullmean_marginal_flat}} → {{nullmean_marginal_conc}} h).
A null that does not track the marginal cannot distinguish shared pile-up from shared
ortholog-specific timing, which is the whole of the effect being measured.

**Calibration under the designs actually used.** Repeating the zero-conservation construction under
the sampling design of each analysed dataset. Across both panels below: {{dm_cells}} cells,
{{dm_outer}} outer trials and {{dm_inner}} inner permutations per cell, {{dm_n}} orthologs,
rejection rate at nominal α = 0.05 with Wilson 95% intervals in parentheses. Phases are drawn
independently for the two datasets, so the true conserved fraction is zero by construction and
every rejection is a false positive; concentration is the fraction of genes whose true phase falls
in one 6 h window. The null hypothesis tested in each cell is that the two phase vectors are
unrelated given the design, and the test is the one-sided empirical tail probability of the
observed mean |phase difference| under the stated null:

{{table:design_table_ours}}

The primary null holds its nominal size under every design this manuscript relies on.

Under the three sampling designs reported by Ferrari et al., with the time-point-permutation null at
the published 1000 permutations:

{{table:design_table_published}}

**This is design-matched, not data-matched.** Phase concentration is swept rather than taken from
their data, so these rows describe how the procedure behaves under conditions resembling their
experiments. They are a property of a sampling design and a null, not of anybody's result. No rate
reported here is attributable to the published analysis, which used different organisms, different
data and a different rhythm caller; attaching any of these numbers to it would be unsupported.

**The coarse design is not monotonic in phase concentration, and must not be reduced to one
inflation factor.** The {{d7_pts}}-point, 4 h design (row `{{d7_design}}` above) gives a
false-positive rate of
{{d7_fpr_c0}} ({{d7_ci_c0}}) at uniform phases, {{d7_fpr_c40}} ({{d7_ci_c40}}) at 40%
concentration, {{d7_fpr_c60}} ({{d7_ci_c60}}) at 60%, and {{d7_fpr_c90}} ({{d7_ci_c90}}) at 90%.
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
{{nullmean_timepoint_flat}} h at uniform phases to {{nullmean_timepoint_conc}} h at 90%
concentration, in the three-null table above. **We have not measured the observed and null locations separately within the
{{d7_pts}}-point design**, so the mechanism is an explanation consistent with the numbers rather
than one this simulation establishes; separating them would require recording both locations per
cell, which these runs do not.

**This also corrects v0.1, which described the published comparison's null as expecting uniform phases.**
It does not. A time-point-permutation null is substantially protected relative to a uniform null,
and it is the appropriate object of criticism. It is nonetheless inflated — a false-positive rate of
{{fpr_timepoint_min}}–{{fpr_timepoint_max}} against a nominal 0.05, and like the uniform null it worsens with more
orthologs ({{fpr_tp_c40_nlo}}→{{fpr_tp_c40_nhi}} at 40% concentration,
{{fpr_tp_c90_nlo}}→{{fpr_tp_c90_nhi}} at 90%, from n={{tn_n_lo}} to n={{tn_n_hi}}) — across the entire range of phase concentration observed in real
data ({{pileup_min}}–{{pileup_max}} here).

### 3.7 Structure, marginals, and the limits of the panel

Genomically adjacent same-strand genes within 2 kb share phase (ratio of adjacent to random mean |Δphase| {{adj_ratio_min}}–{{adj_ratio_max}}; {{adj_within1h_min}}–{{adj_within1h_max}}% of adjacent pairs within 1 h), so genes are not exchangeable
genome-wide. Repeating the tests under circular rotation and operon-block permutation changes
effects and p values negligibly ({{block_gene}} → {{block_rot}} → {{block_blk}}), because the
single-copy ortholog subset is sparse along the genome (§2.5).

Each lineage concentrates transcripts at its own time of day, from {{pileup_min_pct}}% to {{pileup_max_pct}}% in the busiest 6 h
window, with non-coincident windows (Figure 3).

Only two lineages have genome-wide free-running data. The *Halothece* series gives {{halo_timepoints}} time points
over 30 h per condition, {{halo_rhy_min}}–{{halo_rhy_max}}% rhythmic genes, and a best-fit period
pinned at the {{halo_tau}} h scan boundary — underpowered. Effectively one cyanobacterial lineage has adequate constant-light data.

### 3.8 Conditioning on functional composition

Permuting ortholog identity only within broad functional classes asks whether correspondence is
carried by the orthologs themselves or by the functional composition of the sets being compared.

At BRITE level A (6–7 classes):

{{table:strat_table_A}}

At BRITE level B (22–28 classes):

{{table:strat_table_B}}

Three things follow. First, the within-strain benchmark is essentially unaffected at every level
({{stratA_within}}%, {{stratB_within}}%, {{stratB_within_big}}% with small classes dropped), so the
reference scale is not an artifact of functional composition. Second, at coarse classes
{{stratA_sig}} of {{stratA_n}} cross-lineage pairs remain significant, and at fine classes
{{stratB_sig}} of {{stratB_n}} do; with small strata removed only {{stratB_sig_big}} of
{{stratB_n}} do. Third, which pairs those are is not stable, so the aggregate count is the only
level at which we read this analysis; the next paragraph gives the reason and the evidence.

The whole stratified comparison is **exploratory**. It conditions on an annotation hierarchy built
for a different purpose, the pairs are not independent, and the quantity it reports is a count over
seven dependent comparisons, not a test with a controlled error rate.

**Per-pair stratified p values depend on an arbitrary annotation choice.** BRITE is not a partition:
{{multi_pct_A}}% of KOs sit under more than one level-A category and {{multi_pct_B}}% under more
than one level-B category, up to {{multi_max_B}} classes for a single KO. Stratification needs one
class per ortholog, so a rule is required, and the reported analysis takes each KO's first
occurrence in the hierarchy file. Repeating everything with the last occurrence instead — an
equally defensible and equally arbitrary rule — leaves the count of significant cross-lineage pairs
almost unchanged ({{stratA_sig}} → {{stratA_sig_last}} at level A, {{stratB_sig}} →
{{stratB_sig_last}} at level B) but moves {{stratA_n_flip}} pair across α = 0.05 at level A and
{{stratB_n_flip}} at level B. The pair that moves at both levels is *S. elongatus* ×
*Prochlorococcus*: at level A, p = {{stratA_med4_p1}} under the first-occurrence rule and
{{stratA_med4_p2}} under the last, with the effect moving from {{stratA_med4_e1}}% to
{{stratA_med4_e2}}%; at level B, {{stratB_med4_p1}} and {{stratB_med4_p2}}, effect
{{stratB_med4_e1}}% to {{stratB_med4_e2}}%. The pairs that move are, at level A,
{{stratA_flip_pairs}}; at level B, {{stratB_flip_pairs}}. Its effect is small enough that
bookkeeping over multiply-annotated genes decides its significance. We therefore report the
aggregate behaviour of the stratified null and **do not** interpret any individual pair's stratified
p value, including the *Prochlorococcus* one. This supersedes a claim in an earlier draft that the
pair failed under every stratified variant.

What "stable under the assignment rule" means here is narrow and worth stating exactly, because the
word can be stretched to cover three different things. It means: the **count** of cross-lineage
pairs significant at α = 0.05 changes by at most one when the rule is swapped, and the **sign** of
every pair's effect is unchanged. It does **not** mean the identity of the significant pairs is
stable — {{stratA_n_flip}} and {{stratB_n_flip}} pairs cross α — nor that effect **sizes** are
stable, which for small-effect pairs they are not. The first-occurrence rule was fixed as primary
before any stratified p value was examined, and was not chosen for the result it gives; under the
alternative rule the aggregate count would support the same reading.

**The attenuation at fine strata is confounded with power.** Dropping classes smaller than ten
removes {{stratB_drop_lo}}–{{stratB_drop_hi}} orthologs per pair and in two cases roughly halves the
comparison ({{stratB_cut1_pair}} {{stratB_cut1_before}} → {{stratB_cut1_after}};
{{stratB_cut2_pair}} {{stratB_cut2_before}} → {{stratB_cut2_after}}), so loss of significance cannot
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

Aligning lineages by best global shift and ranking orthologs by phase dispersion yields {{core_n}}
tightly co-phased orthologs — dispersion below {{core_cutoff}} h, rhythmic in at least
{{core_min_rhy}} of the five lineages — with biosynthetic annotation (*fabI*, *fabZ*, *dxr*,
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

**Priority.** A documented search ({{lit_queries}} PubMed queries, {{lit_unique}} unique records
screened, {{lit_flagged}} adjudicated; all {{lit_citing}} records citing Ferrari et al. screened at
title and abstract, of which {{lit_citing_cyano}} mention a cyanobacterial taxon and
{{lit_citing_null}} mention a null model) found no study meeting the primary criterion. The strategy
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
3. BRITE assigns {{multi_pct_B}}% of KOs to several classes, and which class a stratification uses
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
- **Per-pair stratified p values withdrawn.** BRITE is not a partition ({{multi_pct_B}}% of KOs
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
