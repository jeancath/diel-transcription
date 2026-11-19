# Prior-art search for Q1

**Question searched:** has anyone already compared diel/circadian transcriptional programs across
cyanobacterial lineages at the level of orthologous genes — which genes cycle, and at what phase?

**Searched:** PubMed (eutils, 11 query formulations across comparative/meta-analysis/conservation/
phase/ortholog axes), plus web search. Abstracts and, for the closest match, full text were read.
Date: 18 Nov 2026.

**Verdict: the comparative design is NOT novel — it was published for plants and algae in 2019, with
the opposite conclusion. The within-cyanobacteria application, the null model, and the within-strain
reference scale do appear to be new. The finding has to be repositioned accordingly.**

---

## 1. The direct methodological precedent — and it concludes the opposite

**Ferrari C, Proost S, Janowski M, et al. "Kingdom-wide comparison reveals the evolution of diurnal
gene expression in Archaeplastida." *Nat Commun* 10:737 (2019). PMID 30760717, PMCID PMC6374488.**

This is essentially the same analysis in a different clade:

| | Ferrari et al. 2019 | this work |
|---|---|---|
| taxa | 9 Archaeplastida (glaucophyte → eudicot), **plus *Synechocystis* PCC 6803 as outgroup** | 6 cyanobacterial lineages |
| orthology | OrthoFinder, one-to-one only | KEGG KO single-copy, **and** Smith-Waterman reciprocal best hits |
| rhythm calling | JTK_Cycle, phase = peak time | cosinor with per-series intercepts, phase = fitted peak |
| statistic | mean \|Δphase\| between orthologs, after scanning a global shift of −12…+12 h | the same |
| null | **permute the time points of the expression matrix, re-run JTK** → near-uniform expected phase | the same, **plus** a null that preserves each dataset's phase distribution and shuffles ortholog identity |
| reference scale | none | within-strain, two independent studies of PCC 7942 |
| conclusion | "diurnal transcriptional programs of these organisms are similar… evidence for the universality of diurnal gene expression" | a weak conserved component, 6–15× smaller than within-strain reproducibility |

Their headline numbers: Δphase between *C. reinhardtii* and *C. paradoxa* peaked at 2–3 h; the pair
was significantly similar at FDR-corrected empirical p < 0.012 from 1000 permutations; most species
pairs became significantly similar once one species' phases were shifted (required shifts up to 9 h);
rhythmic fractions ranged 5.2% (spruce) to 55.9%.

**The caveat worth flagging.** Ferrari et al. report that rhythmic genes pile up — "the majority of
rhythmic genes tend to peak during the second half of the day, at either dusk or night," with
angiosperms bimodal at the two transitions. Their null permutes time points and so expects
near-uniform phases. A statistic compared against a uniform expectation will register "similar"
whenever two organisms share a pile-up, whether or not the timing is ortholog-specific. Preserving
the marginal phase distribution and shuffling only *which ortholog carries which phase* separates
those two explanations.

**Important: tested, not assumed.** Both nulls were run on the cyanobacterial panel
(`src/null_comparison.py`, `results/n1_two_nulls.csv`). On this data they give broadly similar
verdicts (5/10 vs 7/10 cross-lineage pairs significant), so the null choice does **not** by itself
overturn anything in this panel. What separates within-strain from cross-lineage here is effect
size, not the p value.

**But the failure mode is real and quantified.** In a simulation where the truth is fixed at zero
ortholog-specific conservation and each species' phases are drawn independently from its own
piled-up marginal distribution, a uniform-expectation null declares the two transcriptomes
significantly similar **83% of the time when 40% of genes peak within one 6 h window, and 100% of
the time at 60% or 90%** (at 1500 orthologs, 100% from 40% pile-up upward — more orthologs make it
worse, because the added power detects the shared marginal rather than shared ortholog timing);
the marginal-preserving null holds its nominal 5% throughout
(`src/null_simulation.py`). The busiest-6 h fraction in the real datasets analysed here is
0.32–0.91 — all at or above the regime where the uniform null fails at least 83% of the time.

This still does not prove Ferrari et al.'s result is an artifact: their null is a time-point
permutation re-run through JTK, which approaches uniformity without necessarily reaching it (in the
emulation here it gave nulls of 4.4–5.7 h mean |Δphase| against 6.0 h for a uniform distribution).
What it does establish is that the burden of proof now sits with the uniform-style null, and that
re-running the Archaeplastida comparison under a marginal-preserving null is a well-defined and
worthwhile piece of work.

## 2. The mammalian analogue of the within-strain result

**Brooks TG, Manjrekar A, Mrčela A, Grant GR. "Meta-analysis of Diurnal Transcriptomics in Mouse
Liver Reveals Low Repeatability of Rhythm Analyses." *J Biol Rhythms* 38(6):556–570 (2023).
PMID 37382061, PMCID PMC10615793.**

57 mouse-liver time series, 1096 RNA-seq samples. Their findings map almost one-to-one onto the
within-strain analysis here:

* "Core clock genes were remarkably consistent in phase across all studies" → *kaiBC* phase agrees
  to 1.7 h between the two independent 7942 studies.
* "Overlap of genes identified as rhythmic across studies was generally low, with no pair of studies
  having over 60% overlap" → Jaccard 0.549, replication rate 0.690 here.
* "Distributions of phases of significant genes were remarkably inconsistent across studies."

**And a methodological warning that applies directly to this work:** "Despite the discrepancies
between single-study analyses, cross-study analyses found substantial consistency. Running
compareRhythms on each pair of studies identified a median of only 11% of the identified rhythmic
genes as rhythmic in only 1 of the 2 studies." In other words, **thresholding each study separately
manufactures apparent non-replication**. That is precisely the trap the first pass of this analysis
fell into: requiring a gene to clear a significance threshold in *both* datasets, on ~17–175
orthologs per pair, was underpowered and produced a "cross-lineage conservation is at chance"
conclusion that a better-powered continuous analysis then corrected. The correction is recorded in
`REPORT.md`.

## 3. Supporting observations, not the same analysis

* **Kennedy SJ, Risser DD, Paul BG. "Diel expression dynamics in filamentous cyanobacteria."
  *mBio* 16(12):e0377924 (2025). PMID 41251487.** The source of the *Nostoc punctiforme* data used
  here. They report a putative circadian regulator strongly oscillating day/night, "in sharp
  contrast to the arrhythmic pattern observed for a homolog of the canonical circadian input kinase
  in unicellular cyanobacteria" — an independent, single-gene observation of exactly the kind of
  lineage-specific output rewiring this analysis measures genome-wide. No cross-species ortholog
  phase comparison is performed.
* **Hörnlein C, Confurius-Guns V, Stal LJ, Bolhuis H. "Daily rhythmicity in coastal microbial mats."
  *npj Biofilms Microbiomes* 4:11 (2018). PMID 29796291.** Community metatranscriptome, 6 timepoints;
  ~7% of 24,035 conserved gene transcript clusters rhythmic, spanning cyanobacteria, Proteobacteria
  and Bacteroidetes. A community-level question, not an ortholog-level cross-lineage comparison.
* **Schmelling et al., *BMC Evol Biol* 17:169 (2017)** establishes the distribution of clock
  *components* across 11,264 genomes — presence/absence, not output program. The gap this analysis
  addresses is still open in the literature.
* Rhythm-detection methodology: Hughes et al. (JTK_CYCLE), the five-method evaluation in
  *J Biol Rhythms* (2014, PMID 25238853), and the systematic evaluation in *Brief Bioinform* (2021,
  PMID 32672832) justify the choice of a cosinor/F-test caller and the use of a null-strain control.

## 4. What searching did *not* find

No study comparing diel or circadian transcriptomes **between cyanobacterial lineages** at
ortholog level. The nearest things are Ferrari et al. (cross-kingdom, eukaryote-focused, one
cyanobacterium as outgroup) and single-organism diel papers that compare their results to PCC 7942
narratively rather than quantitatively.

## 5. Repositioning

| claim | status after the search |
|---|---|
| "nobody has compared diel programs across species at ortholog level" | **false** — Ferrari et al. 2019 |
| "nobody has done it within cyanobacteria" | appears true |
| "rhythmic gene-set membership replicates poorly between studies" | **known** — Brooks et al. 2023 (mouse); new for cyanobacteria |
| "using within-strain cross-study reproducibility as the yardstick for cross-lineage conservation" | appears new |
| "a marginal-preserving null is required for this statistic" | appears new as an explicit control; a live caveat against Ferrari's uniform-style null, tested but not demonstrated to change the verdict here |
| "the cross-lineage conserved component is ~6–15× smaller than the within-strain effect" | new, and the actual contribution |

**Honest framing for any write-up:** this is a cyanobacterial replication of a published
cross-kingdom design, with a stricter null and a reference scale the original lacked, and it reaches
a quantitatively weaker conclusion about conservation than the original did. That is a worthwhile
contribution and it is not a discovery of a new phenomenon.

## 6. Search limitations

PubMed title/abstract queries plus one extended web search; no Google Scholar citation-graph
traversal, no screening of the reference lists of Ferrari et al. or the mBio paper, no bioRxiv
full-text search. A study that performed this comparison incidentally, inside a paper about
something else, would likely have been missed. Before submission, the citing literature of
Ferrari et al. 2019 should be screened directly.
