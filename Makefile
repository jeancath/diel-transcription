# Single entry point. Inputs download to $(DIEL_DATA_DIR), default ./data .
#
# `make all` regenerates every file the manuscript reads, from the raw GEO
# series onward. It must stay that way: `make clean` deletes all derived
# output, so anything `all` does not rebuild would be lost. The semantic
# validator enforces this (every results/ file must have a producer in src/).
#
# The Monte Carlo targets dominate the runtime (~1-2 h). They checkpoint per
# grid cell into results/*.partial.csv and resume, so an interrupted run does
# not start over.
PY ?= python3
.PHONY: all data rhythms analysis robustness nulls simulations strat literature \
        orthology figures clean clean-tests check manuscript check-manuscript validate tests-slow

# Order matters: `nulls` reads the genome annotations that `data` fetches, and
# `manuscript` reads outputs of every stage before it. A clean-checkout build
# failed here once because orthology came after nulls.
all: rhythms analysis robustness orthology nulls simulations strat literature \
     figures manuscript check-manuscript validate

data:          ; $(PY) src/fetch_data.py
rhythms: data  ; $(PY) src/run_rhythms.py
analysis:      ; $(PY) src/analyze.py
robustness:    ; $(PY) src/robustness.py
nulls:         ; $(PY) src/null_comparison.py
	$(PY) src/threshold_sensitivity.py
	$(PY) src/independence_check.py
	$(PY) src/block_null.py
	$(PY) src/dependence_of_pairs.py
	$(PY) src/effect_vs_pileup.py
	$(PY) src/conserved_core.py

# Monte Carlo. The four long jobs are independent; run them concurrently and
# fail the target if any of them fails.
simulations:   ; $(PY) src/null_simulation.py
	$(PY) src/power_analysis.py
	$(PY) src/power_grid.py
	@set -e; pids=""; \
	for s in three_nulls_fast power_grid_full power_small_fractions design_matched_nulls; do \
	  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 $(PY) src/$$s.py > results/$$s.log 2>&1 & \
	  pids="$$pids $$!"; done; \
	rc=0; for p in $$pids; do wait $$p || rc=1; done; \
	if [ $$rc -ne 0 ]; then echo "a simulation failed; see results/*.log"; exit 1; fi

# functional-stratified null plus its closure report at both BRITE levels
strat:         ; $(PY) src/stratified_null.py
# prior-art search and the citation screen of the closest prior study
literature:    ; $(PY) src/literature_screen.py
# Smith-Waterman reciprocal best hits. Slowest non-Monte-Carlo stage (~1 h);
# in `all` because `clean` removes its output and the manuscript reads it.
orthology:     ; $(PY) src/rbh.py --subset
	$(PY) src/rbh_check.py
figures:       ; $(PY) src/figures.py

manuscript:    ; $(PY) src/manuscript_numbers.py > /dev/null
	$(PY) src/render_manuscript.py
# fails if any statistic in the manuscript was typed by hand or has drifted
check-manuscript: ; $(PY) src/check_manuscript.py
# semantic assertions: identity, provenance, claim strength, and reproducibility
validate:      ; $(PY) src/validate_semantics.py
	$(PY) src/audit_renderer.py
	$(PY) src/test_cosinor_equivalence.py
	$(PY) src/test_semantic_mutations.py
	$(PY) src/coverage_report.py
	$(PY) src/number_provenance.py
	$(PY) src/provenance.py
# the slow behavioural tests, not in `all`: interrupt/resume determinism
tests-slow:    ; $(PY) src/test_determinism.py
check:         ; $(PY) -c "import numpy,pandas,scipy,matplotlib,Bio,pyswrd; print('deps ok')"
# removes everything derived; leaves data/ (the raw inputs) alone
# `clean` removes exactly what `make all` regenerates, so that clean + all is a
# round trip. The determinism artifacts come from tests-slow, which is not in
# all and takes much longer, so they are cleaned separately and deliberately.
clean:         ; rm -f results/*.csv results/*.csv.gz results/*.log figures/*.png
	find results -maxdepth 1 -name '*.json' ! -name 'determinism_report.json' -delete
clean-tests:   ; rm -rf results/_determinism results/determinism_report.json
