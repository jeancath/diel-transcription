# Literature search strategy and screening record

Supporting the priority statement in `MANUSCRIPT.md`. Everything here is
reproducible: `python3 src/literature_screen.py` re-runs the query set and
regenerates `results/lit_query_log.csv` and `results/lit_screen.csv`.

**Search dates:** 18 November 2026 (all searches).
**Databases:** PubMed/MEDLINE via NCBI E-utilities; PubMed citation network via
`elink` (`pubmed_pubmed_citedin`); one extended general web search.
**Not searched:** bioRxiv/medRxiv full text, Google Scholar citation graph, Web of
Science, Scopus, EMBASE. This is a stated limitation, not an omission.

---

## 1. Review question

**Primary.** Has any study compared diel or circadian transcriptional programs
**between two or more cyanobacterial lineages at the level of orthologous genes**
— that is, which orthologs cycle and at what phase?

**Secondary.** Has any study addressed the **null model** used when comparing
event timing across species, specifically the effect of non-uniform marginal
phase distributions?

## 2. Inclusion and exclusion criteria

| | Include | Exclude |
|---|---|---|
| organism | two or more cyanobacterial lineages | single strain; community metatranscriptome without lineage-level ortholog comparison; non-cyanobacterial only |
| measurement | genome-wide transcriptome time course | single-gene or reporter-only; proteome-only |
| comparison | ortholog-level comparison of which genes cycle and/or their phase | narrative comparison to published work without a joint analysis; comparison of clock-component presence/absence only |
| outcome | rhythmic gene set and/or phase correspondence | growth, yield, or physiology only |

## 3. Searches and yields

### 3.1 Query set (PubMed, title/abstract)

| id | records retrieved | records screened |
|---|---|---|
| Q1 cyanobacteria × comparative × diel transcriptome | 15 | 15 |
| Q2 cyanobacteria × rhythm × meta-analysis/reanalysis | 0 | 0 |
| Q3 clock output × cyanobacteria × conservation/evolution | 9 | 9 |
| Q4 circadian/diel × phase × ortholog × conservation | 14 | 14 |
| Q5 null model/permutation × phase/timing × species/ortholog | 2 | 2 |
| Q6 circadian/diurnal transcriptome × cross-species | 39 | 39 |
| Q7 circadian/rhythmic × reproducibility × transcriptome | 187 | 187 |
| Q8 *Synechococcus elongatus* × circadian × transcriptome | 5 | 5 |
| **unique records after deduplication** | | **265** |

Exact query strings are in `src/literature_screen.py` and reproduced in
`results/lit_query_log.csv`.

### 3.2 Citation network of the closest prior study

Ferrari et al., *Nat Commun* 10:737 (2019), PMID 30760717 — the methodological
precedent — was cited by **34 records** (PubMed `citedin`, 18 Nov 2026).
**All 34 were screened at title and abstract level** (`results/lit_ferrari_citing_screen.csv`).

* records mentioning any cyanobacterial taxon in title or abstract: **0**
* records mentioning a null model, permutation or background distribution: **0**
* records with timing terms: 24; with comparative/ortholog terms: 14

The citing literature is almost entirely plant and algal. The nearest items —
a cross-species *Arabidopsis*/soybean circadian comparison (PMID 37836084) and a
potato diel conservation study (PMID 41102621) — were read at abstract level;
neither compares cyanobacterial lineages and neither addresses the null model.

### 3.3 Screening of the 265 unique records

Automated title/abstract screening against the criteria flagged **31** records
for adjudication. All were adjudicated by reading title and abstract; four were
read in further detail. None met the primary criterion. Representative
adjudications:

| record | why it does not meet the primary criterion |
|---|---|
| Schmelling et al. 2017, PMID 28732467 | clock **component** presence/absence across 11,264 genomes; no expression data |
| Hörnlein et al. 2018, PMID 29796291 | community metatranscriptome of a microbial mat; no lineage-level ortholog phase comparison |
| Muñoz-Marín et al. 2023, PMID 37130101 | two strains of one UCYN-A species; transcriptome comparison, not diel phase correspondence |
| Kennedy et al. 2025, PMID 41251487 | single-organism *Nostoc punctiforme* diel transcriptome; compares to PCC 7942 narratively, not by joint analysis |
| Ferrari et al. 2019, PMID 30760717 | Archaeplastida; one cyanobacterium as outgroup, not compared with other cyanobacteria |
| various | single-organism time courses, clock mechanism, biotechnology applications, non-cyanobacterial clocks |

Two records flagged by the secondary criterion were regex false positives on
"phase" and "permutation" (a spider survey, PMID 23321102; an invasive-plant
climate model, PMID 40817242). They are recorded rather than silently dropped,
because they show the automated screen is loose and that adjudication by reading
is doing real work.

### 3.4 Methodologically relevant work found by the search

Not prior art for the primary question, but it belongs in the manuscript:

* **Lehmann et al. 2013, *BMC Bioinformatics* 14:133, PMID 23601192** — "How
  cyanobacteria pose new problems to old methods: challenges in microarray time
  series analysis." Directly about analytical pitfalls in cyanobacterial diurnal
  transcriptome time series, including that interpretation usually focuses on
  time and phase of expression. This should be cited.
* **Brooks et al. 2023, PMID 37382061** — mouse-liver meta-analysis; the
  within-study/between-study dissociation this work reproduces in a bacterium.
* **Mei et al. 2021, PMID 32672832; Hughes et al. 2014, PMID 25238853** —
  evaluations of genome-wide rhythm-detection methods.

## 4. Conclusion and the wording it licenses

No study meeting the primary criterion was found. The search does **not**
establish that none exists: the citation graph was traversed one level from a
single anchor paper, preprint full text was not searched, and a comparison made
incidentally inside a paper about something else would have been missed.

The manuscript therefore states:

> To our knowledge, this question has not been tested within cyanobacterial
> lineages.

and does not claim priority, novelty, or that the comparison is the first of its
kind. The contribution rests on the comparison framework and the within-strain
benchmark, not on being first.

## 5. Known limitations of this search

1. One-level citation traversal from one anchor paper.
2. No preprint full-text search.
3. PubMed only; no Web of Science, Scopus or Google Scholar.
4. Automated screening by regular expression, with human adjudication only of
   flagged records — a relevant record whose title and abstract avoid all the
   search terms would not have surfaced.
5. English-language records only, implicitly, since no language filter was
   applied but PubMed indexing dominates.
