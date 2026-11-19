# Within-strain reproducibility sets the scale for cross-species comparison of diel transcriptional programs

Melad Naieem Fathy · Independent researcher ·
[ORCID 0009-0002-4553-8894](https://orcid.org/0009-0002-4553-8894)

Code, derived data, results and manuscript for a reanalysis of public cyanobacterial diel
transcriptome series. This is the public release of a private research repository; it is an
allowlisted snapshot, not a mirror of the working history.

<!-- DOI badge: added after the Zenodo release, using the CONCEPT DOI (all versions). -->

## The question

Diel transcriptional programs in cyanobacteria were mapped principally in one strain,
*Synechococcus elongatus* PCC 7942, and are often discussed as if representative of the group.
Does the identity and timing of rhythmic orthologs transfer between cyanobacterial lineages?

## What this study adds

**A benchmark the field was missing.** Cross-species comparisons of transcript timing report how
much better than chance two species agree, but nothing said how much two studies of the *same*
organism agree. Two independent studies of one strain, from different laboratories on different
platforms, supply that scale. Within that strain ortholog peak times agree to a median of 1.67 h.
Across lineages the same statistic is far smaller — a median 5.0% reduction against a
phase-preserving null, against 75.0% within the strain.

**Null calibration.** On simulated data with conservation fixed at zero, a null that expects
uniform phases has a false-positive rate of 0.80–1.00 under realistic phase concentration, and a
time-point-permutation null of the kind used in a published cross-species comparison is partially
but not fully protected (0.21–0.92). Only a permutation of ortholog identity that preserves each
dataset's phase distribution holds near nominal (0.04–0.06). Every dataset examined here lies in
the regime where the first two are inflated.

**A bound, with its conditions stated.** Power analysis shows the design excludes timing
conservation involving 40% of orthologs under tight phase jitter, and that this exclusion is
substantially weaker under 4 h jitter. 40% is not a majority of the program, and the bound is
conditional on the simulation's assumptions. Both limits are stated in the paper rather than left
to the reader.

## What this study does not show

- **It is about diel transcription, not the circadian clock.** The data are from cells under
  light/dark cycles. A diel transcript rhythm can be driven by the light/dark cycle, gated by an
  endogenous oscillator, or both, and these data cannot separate them. Nothing here settles whether
  endogenous circadian regulation is conserved.
- **"Not detected" is not "absent".** Where a comparison shows no signal, that is what is reported.
- **It is exploratory throughout.** No pre-registration was frozen, and several central analysis
  choices were made with the data in view. `prereg/NO-PREREGISTRATION.md` lists which ones and when.
- **The comparative design is not new.** Ferrari et al. 2019 (*Nat Commun* 10:737) ran it for
  Archaeplastida and reached the opposite conclusion. `PRIOR-ART.md` says so before any novelty
  claim, and the null-calibration result is offered as an explanation of how two such studies can
  disagree — not as a verdict on that paper.

## Reading order

| File | What it is |
|---|---|
| `MANUSCRIPT.md` | the paper; every statistic in it is generated from `results/` |
| `CLOSURE-REPORT.md` | three separate verdicts: content reproducibility, byte-level reproducibility, scientific validity |
| `PRIOR-ART.md` | read before citing anything here as new |
| `REPORT.md` | the longer technical report |
| `DEVIATIONS.md` | numbered, dated record of every change of analysis or claim |
| `RELEASE-v0.3.md` | every released file with a SHA-256 |

## Reproducing it

```bash
pip install -r requirements.txt
make data                      # downloads the public GEO/NCBI inputs; none are redistributed here
make clean && make check && make all
python3 audit/number_check.py
```

The primary expression data are not included in this repository — they are public GEO series and
NCBI records, and `src/fetch_data.py` (run by `make data`) downloads every one of them without
credentials. Until it has run, the check that every accession resolves to a downloaded series
matrix will fail; that is expected in a fresh clone, not a defect.

`make check` runs the semantic assertions, the mutation tests, the renderer audit and the
hand-typed-statistics scan, and fails the build on any of them. The full run recomputes all 93
simulation cells; it is single-threaded by design, because each cell's random number generator is
seeded from the cell's coordinates rather than drawn from one advancing stream, which is what makes
the result reproducible.

Three independent full recomputations agree bit-for-bit. Content-level reproducibility is complete;
byte-level reproducibility is short by exactly one file, `results/rhythms_per_gene.csv.gz`, because
gzip writes the current time into its header.

## Use of AI tools

The author found the idea, critiqued it, searched the literature, reviewed the work, and analysed
and interpreted the results. A generative AI assistant also did some of this work, including
writing code.

## Licence

MIT for code; CC BY 4.0 for text, figures and data tables. See `LICENSE`. The primary expression
data are public GEO and NCBI records and remain under their depositors' terms; only derived tables
are redistributed here.

## Citation

See `CITATION.cff`.
