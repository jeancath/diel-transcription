"""Single source of truth for every number quoted in the manuscript.

Each entry is computed from a file in results/ at call time. MANUSCRIPT.md is
rendered from MANUSCRIPT.template.md by substituting {{key}} placeholders with
these values, and src/check_manuscript.py fails if the rendered file and the
template disagree. No headline number is typed by hand.
"""
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
WITHIN = {"SE7942_ito_WT", "SE7942_vijayan"}
LAB = {"SE7942_ito_WT": "*S. elongatus*", "SE7942_vijayan": "*S. elongatus* (study 2)",
       "SYN6803": "*Synechocystis*", "TERY": "*Trichodesmium*",
       "MED4": "*Prochlorococcus*", "NPUN": "*Nostoc*"}


# which make target produces which result file, so a missing input names its
# own fix instead of a bare traceback
_PRODUCER = {
    "rhythms_per_gene": "rhythms", "dataset_summary": "rhythms",
    "q1_": "analysis", "q2_": "analysis", "q3_": "analysis",
    "clock_gene_control": "analysis",
    "r1_": "robustness", "r2_": "robustness", "r3_": "robustness",
    "n1_": "nulls", "n2_": "nulls", "n5_": "nulls", "n7_": "nulls",
    "n8_": "nulls", "n10_": "nulls", "n11_": "nulls",
    "n3_": "orthology",
    "n4_": "simulations", "n6_": "simulations",
    "n12_": "simulations", "n13_": "simulations", "n14_": "simulations",
    "n16_": "simulations", "n17_": "simulations",
    "n15_": "strat", "n18_": "strat", "n19_": "strat",
    "lit_": "literature",
    "determinism_report": "tests-slow",
    "numbers": "manuscript", "fig1_": "figures", "fig2_": "figures",
    "fig3_": "figures",
    "PROVENANCE": "validate", "ASSERTION_COVERAGE": "validate",
    "NUMBER_PROVENANCE": "validate",
}


def producer(name):
    """The make target that writes this result file, or None if unknown."""
    return next((t for k, t in _PRODUCER.items() if name.startswith(k)), None)


def _csv(name):
    f = OUT / name
    if not f.exists():
        tgt = producer(name)
        hint = f"run `make {tgt}`" if tgt else "no known producer - check the Makefile"
        raise FileNotFoundError(f"{name} is missing from results/; {hint}")
    return pd.read_csv(f)


def collect():
    v = {}

    # ---- panel and controls
    ds = _csv("dataset_summary.csv")
    v["n_datasets"] = len(ds)            # analysed arms, incl. the two Halothece salinities
    v["n_series"] = int(ds.gse.nunique())  # distinct GEO series
    v["n_lineages"] = ds.lineage.str.replace(r" \(.*\)| 0\.5M| 2\.5M", "", regex=True).nunique()
    kn = ds[ds.dataset == "SE7942_ito_KAINULL"].iloc[0]
    v["kainull_rhythmic"] = int(kn.n_rhythmic)
    v["kainull_genes"] = int(kn.n_genes)
    cg = _csv("clock_gene_control.csv")
    for g in ("kaiB", "kaiC"):
        k = cg[(cg.clock_gene == g) & (cg.dataset == "SE7942_ito_WT")].iloc[0]
        v[f"{g}_amp"] = round(float(k.amp_log2), 2)
        # the locus and four significant figures, so a reader can see that two
        # equal-looking amplitudes come from two different genes
        v[f"{g}_amp4"] = round(float(k.amp_log2), 4)
        v[f"{g}_locus"] = str(k.gene)
        v[f"{g}_phase"] = round(float(k.phase_h), 1)
        v[f"{g}_q"] = f"{float(k.q):.1e}"
    kt = cg[(cg.clock_gene == "kaiC") & (cg.lineage == "Trichodesmium")].iloc[0]
    v["kaiC_tery_phase"] = round(float(kt.phase_h), 1)

    # ---- within-strain reproducibility
    q1 = _csv("q1_within_species_reproducibility.csv").set_index("metric")["value"]
    v["within_genes"] = int(q1["genes compared (both studies)"])
    v["within_rhy_1"] = int(q1["rhythmic, Ito 2009 (LL, GPL8012)"])
    v["within_rhy_2"] = int(q1["rhythmic, Vijayan 2009 (LL, GPL9534)"])
    v["within_rhy_both"] = int(q1["rhythmic in both"])
    v["within_rhy_exp"] = round(float(q1["expected overlap if independent"]), 0).__int__()
    v["within_median_dphase"] = round(float(q1["median |phase difference| (h)"]), 2)
    v["within_pct_4h"] = f"{100*float(q1['phase agreement within 4 h']):.1f}"
    v["within_jaccard"] = round(float(q1["Jaccard of rhythmic sets"]), 3)

    # ---- effect sizes by selection rule
    ts = _csv("n2_threshold_sensitivity.csv")
    ts["is_within"] = ts.apply(lambda r: {r.a, r.b} == WITHIN, axis=1)
    for rule, tag in (("all orthologs", "all"), ("either q<0.05", "either"),
                      ("both q<0.05", "both"), ("top 200 combined", "top200")):
        g = ts[ts.rule == rule]
        w = g[g.is_within]
        c = g[~g.is_within]
        v[f"within_{tag}"] = round(float(w.effect_pct.iloc[0]), 1)
        v[f"cross_{tag}_median"] = round(float(c.effect_pct.median()), 1)
        v[f"cross_{tag}_max"] = round(float(c.effect_pct.max()), 1)
        v[f"cross_{tag}_sig"] = int((c.p < .05).sum())
    b = ts[(ts.rule == "both q<0.05")]
    w = b[b.is_within].iloc[0]
    v["within_abs_gain_h"] = round(float(w.null_mean - w.obs_mean_dphase_h), 2)

    # ---- leave-one-out / pair dependence
    loo = _csv("n10_leave_one_out.csv")
    r = loo[loo.rule == "both q<0.05"].iloc[0]
    v["loo_all"] = round(float(r.all_pairs), 1)
    for col, tag in (("drop_S.", "drop_SE"), ("drop_Nostoc", "drop_NPUN"),
                     ("drop_Prochlorococcus", "drop_MED4")):
        v[tag] = round(float(r[col]), 1)
    pv = _csv("n11_effect_vs_pileup.csv")
    pv = pv[(~pv.apply(lambda x: {x.a, x.b} == WITHIN, axis=1)) & (pv.rule == "both q<0.05")]
    pv = pv.assign(gain=pv.null_mean - pv.obs_mean_dphase_h,
                   se=pv.a.str.startswith("SE7942") | pv.b.str.startswith("SE7942"))
    v["gain_with_SE_h"] = round(float(pv[pv.se].gain.median()), 2)
    v["gain_without_SE_h"] = round(float(pv[~pv.se].gain.median()), 2)

    # ---- effect vs pile-up correlations (exploratory; pairs are not independent)
    from scipy import stats as _st
    ep = _csv("n11_effect_vs_pileup.csv")
    ep = ep[~ep.apply(lambda r: {r.a, r.b} == WITHIN, axis=1)]
    rp, rn = [], []
    for _rule, _g in ep.groupby("rule"):
        rp.append(_st.spearmanr(_g.pile_mean, _g.effect_pct)[0])
        rn.append(_st.spearmanr(_g.null_mean, _g.effect_pct)[0])
    v["pileup_rho"] = f"{min(rp):+.2f} to {max(rp):+.2f}"
    v["nullmean_rho"] = f"{min(rn):+.2f} to {max(rn):+.2f}"

    # ---- marginals
    pm = _csv("r3_phase_marginals.csv")
    pm = pm[pm.n > 500]
    v["pileup_min"] = f"{pm.frac_in_busiest_6h.min():.2f}"
    v["pileup_max"] = f"{pm.frac_in_busiest_6h.max():.2f}"
    v["pileup_min_pct"] = f"{100*pm.frac_in_busiest_6h.min():.0f}"
    v["pileup_max_pct"] = f"{100*pm.frac_in_busiest_6h.max():.0f}"

    # ---- independence / block nulls
    ic = _csv("n7_independence_check.csv")
    v["adj_ratio_min"] = round(float(ic.ratio.min()), 2)
    v["adj_ratio_max"] = round(float(ic.ratio.max()), 2)
    v["adj_within1h_min"] = f"{100*ic.frac_adjacent_within_1h.min():.0f}"
    v["adj_within1h_max"] = f"{100*ic.frac_adjacent_within_1h.max():.0f}"
    bn = _csv("n8_block_null.csv")
    row = bn[bn.b == "SYN6803"].iloc[0]
    v["block_gene"] = f"{row.effect_gene:.1f}%/p={row.p_gene:.3f}"
    v["block_rot"] = f"{row.effect_rotation:.1f}%/p={row.p_rotation:.3f}"
    v["block_blk"] = f"{row.effect_block:.1f}%/p={row.p_block:.3f}"
    v["mean_block_size"] = round(float((bn.n / bn.n_blocks).mean()), 1)

    # ---- RBH
    rb = _csv("n3_rbh_orthologs.csv")
    rb = rb[~rb.apply(lambda r: {r.a, r.b} == WITHIN, axis=1)]   # drop the identity case
    v["rbh_min"] = int(rb.n_rbh.min())
    v["rbh_max"] = int(rb.n_rbh.max())

    # ---- three nulls
    #
    # There was a fallback here to n9_three_nulls.csv, the 50-trial pilot. It
    # could never fire in a clean build: three_nulls.py is in no make target,
    # so `make clean` deleted n9 and `make all` did not restore it. The branch
    # existed only to make a stale tracked file look like a live input, and it
    # would have quietly published 50-trial numbers if n13 were ever missing.
    # n13 is now required and a missing file is an error.
    tn = _csv("n13_three_nulls_1000.csv")
    v["tn_outer"] = int(tn.outer_trials.iloc[0])
    v["tn_inner"] = int(tn.inner_permutations.iloc[0])
    cols = {"uniform": "FPR_uniform", "timepoint": "FPR_timepoint",
            "marginal": "FPR_marginal"}
    v["tn_source"] = "n13_three_nulls_1000.csv"
    sub = tn[tn.frac_in_6h_window > 0]
    for tag, col in cols.items():
        v[f"fpr_{tag}_min"] = f"{sub[col].min():.2f}"
        v[f"fpr_{tag}_max"] = f"{sub[col].max():.2f}"
    # the mechanism: how far each null's own mean |dphase| moves as phase
    # concentration rises. Taken from the largest n, flat to concentrated.
    if "null_mean_uniform_h" in tn.columns:
        big = tn[tn.n_orthologs == tn.n_orthologs.max()].sort_values("frac_in_6h_window")
        for tag in ("uniform", "timepoint", "marginal"):
            c = f"null_mean_{tag}_h"
            v[f"nullmean_{tag}_flat"] = f"{big[c].iloc[0]:.2f}"
            v[f"nullmean_{tag}_conc"] = f"{big[c].iloc[-1]:.2f}"
    v["tn_cells"] = int(len(tn))
    v["pileup_sim_max"] = f"{tn.frac_in_6h_window.max():.0%}"

    # ---- literature screen
    try:
        lq = _csv("lit_query_log.csv")
        ls = _csv("lit_screen.csv")
        lf = _csv("lit_ferrari_citing_screen.csv")
        v["lit_queries"] = len(lq)
        v["lit_unique"] = len(ls)
        v["lit_flagged"] = int(ls.meets_primary.sum())
        v["lit_citing"] = len(lf)
        v["lit_citing_cyano"] = int(lf.cyano.sum())
        v["lit_citing_null"] = int(lf.nullmodel.sum())
    except Exception:
        pass

    # ---- functional stratification headline
    f = OUT / "n18_stratified_closure_A.csv"
    if f.exists():
        sa = pd.read_csv(f)
        cross = sa[~sa.apply(lambda r: {r.a, r.b} == WITHIN, axis=1)]
        v["stratA_sig"] = int((cross.all_strat_p < .05).sum())
        v["stratA_n"] = len(cross)
        w = sa[sa.apply(lambda r: {r.a, r.b} == WITHIN, axis=1)].iloc[0]
        v["stratA_within"] = round(float(w.all_strat_effect), 1)
    f = OUT / "n18_stratified_closure_B.csv"
    if f.exists():
        sb = pd.read_csv(f)
        cross = sb[~sb.apply(lambda r: {r.a, r.b} == WITHIN, axis=1)]
        v["stratB_sig"] = int((cross.all_strat_p < .05).sum())
        v["stratB_sig_big"] = int((cross.big_classes_only_strat_p < .05).sum())
        v["stratB_n"] = len(cross)
        w = sb[sb.apply(lambda r: {r.a, r.b} == WITHIN, axis=1)].iloc[0]
        v["stratB_within"] = round(float(w.all_strat_effect), 1)
        v["stratB_within_big"] = round(float(w.big_classes_only_strat_effect), 1)
        # how much of each comparison the small-stratum filter removes, and the
        # two pairs it cuts hardest - the confound with power, stated in numbers
        v["stratB_drop_lo"] = int(cross.n_in_classes_lt10.min())
        v["stratB_drop_hi"] = int(cross.n_in_classes_lt10.max())
        worst = cross.assign(kept=cross.big_classes_only_n / cross.all_n) \
                     .sort_values("kept").head(2)
        for i, (_, r) in enumerate(worst.iterrows(), start=1):
            v[f"stratB_cut{i}_before"] = int(r.all_n)
            v[f"stratB_cut{i}_after"] = int(r.big_classes_only_n)
            v[f"stratB_cut{i}_pair"] = f"{LAB.get(r.a, r.a)} \u00d7 {LAB.get(r.b, r.b)}"

    # ---- BRITE is not a partition, and the rule for multiply-assigned KOs
    #      changes per-pair verdicts. Reported, not buried.
    ma = OUT / "n19_brite_multiassign.csv"
    if ma.exists():
        m = pd.read_csv(ma).set_index("brite_level")
        for lvl in ("A", "B"):
            if lvl in m.index:
                v[f"multi_pct_{lvl}"] = m.loc[lvl, "pct_multiassigned"]
        v["multi_max_B"] = int(m.loc["B", "max_classes_per_ko"])
    for lvl in ("A", "B"):
        f1 = OUT / f"n18_stratified_closure_{lvl}.csv"
        f2 = OUT / f"n18_stratified_closure_{lvl}_last.csv"
        if f1.exists() and f2.exists():
            t1, t2 = pd.read_csv(f1), pd.read_csv(f2)
            mg = t1.merge(t2, on=["a", "b"], suffixes=("_1", "_2"))
            cr = mg[~mg.apply(lambda r: {r.a, r.b} == WITHIN, axis=1)]
            v[f"strat{lvl}_sig_last"] = int((cr.all_strat_p_2 < .05).sum())
            flip = cr[(cr.all_strat_p_1 < .05) != (cr.all_strat_p_2 < .05)]
            v[f"strat{lvl}_n_flip"] = int(len(flip))
            v[f"strat{lvl}_flip_pairs"] = ("; ".join(
                f"{LAB.get(r.a, r.a)} \u00d7 {LAB.get(r.b, r.b)} "
                f"(p {r.all_strat_p_1:.3f} vs {r.all_strat_p_2:.3f})"
                for _, r in flip.iterrows()) or "none")
            med = cr[(cr.b == "MED4") | (cr.a == "MED4")]
            if len(med):
                r = med.iloc[0]
                v[f"strat{lvl}_med4_p1"] = f"{r.all_strat_p_1:.3f}"
                v[f"strat{lvl}_med4_p2"] = f"{r.all_strat_p_2:.3f}"
                v[f"strat{lvl}_med4_e1"] = f"{r.all_strat_effect_1:.1f}"
                v[f"strat{lvl}_med4_e2"] = f"{r.all_strat_effect_2:.1f}"

    # ---- power
    #
    # Read from the full surface (n14 + n17), not from n12_power_grid.csv.
    # n12 is a 40-trial pilot; quoting it in the prose while the table beside
    # it reported the 400-trial surface meant two different experiments were
    # being described as one. The trial counts come from the files.
    pg = pd.concat([_csv("n14_power_surface.csv"), _csv("n17_power_small_fractions.csv")],
                   ignore_index=True)
    v["power_outer"] = int(pg.outer_trials.iloc[0])
    v["power_inner"] = int(pg.inner_permutations.iloc[0])
    v["power_cells"] = int(len(pg))
    tight = pg[(pg.true_conserved_frac == 0.4) & (pg.jitter_sd_h <= 2)]
    loose = pg[(pg.true_conserved_frac == 0.4) & (pg.jitter_sd_h == 4)]
    v["power40_tight_min"] = f"{tight.power.min():.2f}"
    v["power40_tight_max"] = f"{tight.power.max():.2f}"
    v["power40_loose_min"] = f"{loose.power.min():.2f}"
    v["power40_loose_max"] = f"{loose.power.max():.2f}"
    small = pg[pg.true_conserved_frac.isin([0.1, 0.2])]
    v["power_small_min"] = f"{small.power.min():.2f}"
    v["power_small_max"] = f"{small.power.max():.2f}"
    null = pg[pg.true_conserved_frac == 0.0]
    v["power_size_min"] = f"{null.power.min():.2f}"
    v["power_size_max"] = f"{null.power.max():.2f}"

    # ---- cross-lineage effect excluding Prochlorococcus (MED4)
    #
    # The quoted range is the range of the per-rule MEDIAN across the four
    # selection rules, not the range of individual pair effects. It was a
    # hand-typed literal until the numeric-literal check was repaired.
    pb = _csv("n10_pairs_bh_corrected.csv")
    pb = pb[pb.a != "a"].copy()
    for c in ("effect_pct", "q_bh"):
        pb[c] = pb[c].astype(float)
    xl = pb[pb.comparison == "cross-lineage"]
    nm = xl[~(xl.a.eq("MED4") | xl.b.eq("MED4"))]
    meds = nm.groupby("rule").effect_pct.median()
    v["exmed4_median_min"] = f"{meds.min():.1f}"
    v["exmed4_median_max"] = f"{meds.max():.1f}"
    both = nm[nm.rule == "both q<0.05"]
    v["exmed4_nsig"] = int((both.q_bh < 0.05).sum())
    v["exmed4_npairs"] = int(len(both))

    # ---- time-point-permutation null: growth of the false-positive rate with n
    tn3 = _csv("n13_three_nulls_1000.csv")
    lo, hi = int(tn3.n_orthologs.min()), int(tn3.n_orthologs.max())
    v["tn_n_lo"], v["tn_n_hi"] = lo, hi

    def _tp(frac, n):
        r = tn3[(tn3.frac_in_6h_window == frac) & (tn3.n_orthologs == n)]
        return f"{float(r.FPR_timepoint.iloc[0]):.2f}"
    for frac, tag in ((0.4, "c40"), (0.9, "c90")):
        v[f"fpr_tp_{tag}_nlo"] = _tp(frac, lo)
        v[f"fpr_tp_{tag}_nhi"] = _tp(frac, hi)

    # ---- the one free-running series other than PCC 7942
    ds = _csv("dataset_summary.csv")
    halo = ds[ds.dataset.str.startswith("HALO")]
    v["halo_rhy_min"] = f"{100 * halo.frac_rhythmic.min():.1f}"
    v["halo_rhy_max"] = f"{100 * halo.frac_rhythmic.max():.1f}"
    v["halo_timepoints"] = int(halo.n_timepoints.max())
    v["halo_tau"] = f"{float(halo.tau_fit.max()):.0f}"

    # ---- the coarse published design, whose false-positive rate is NOT
    # monotonic in phase concentration. Reported cell by cell with intervals,
    # never summarised as one inflation factor.
    dm = _csv("n16_design_matched_nulls.csv")
    v["dm_cells"] = int(len(dm))
    v["dm_outer"] = int(dm.outer_trials.max())
    v["dm_inner"] = int(dm.inner_permutations.max())
    v["dm_n"] = int(dm.n_orthologs.max())
    d7 = dm[dm.design.str.startswith("F: 7 pts")]
    if len(d7):
        v["d7_design"] = str(d7.design.iloc[0])
        v["d7_pts"] = int(d7.n_timepoints.iloc[0])
        for frac in (0.0, 0.4, 0.6, 0.9):
            r = d7[d7.frac_in_6h_window == frac]
            if len(r):
                tag = f"c{int(frac * 100)}"
                v[f"d7_fpr_{tag}"] = f"{float(r.FPR.iloc[0]):.3f}"
                v[f"d7_ci_{tag}"] = str(r.FPR_ci95.iloc[0])

    # ---- the candidate phase-conserved set (3.9)
    #
    # Only the descriptive numbers are published: the set size, the dispersion
    # cutoff and the membership criterion. The null mean, fold enrichment and
    # p value from the fixed-alignment identity shuffle are deliberately NOT
    # exposed as registry keys and so cannot appear in the manuscript. The
    # alignment shift is fitted on the observed ortholog pairings and held
    # fixed across permutations, so that null does not reproduce the selection
    # procedure and cannot calibrate this statistic. The full diagnostic stays
    # in results/n5_core_sd_null.csv, which the manuscript cites as a record
    # rather than as a test. If a future version refits the shift inside each
    # permutation, that is the point at which a p value may be published.
    sn = _csv("n5_core_sd_null.csv")
    m = sn[sn.min_rhythmic == sn.min_rhythmic.max()].iloc[0]
    v["core_n"] = int(m.observed_below_cutoff)
    v["core_cutoff"] = f"{float(m.sd_cutoff_h):.2f}"
    v["core_min_rhy"] = int(m.min_rhythmic)
    return v


def tables(vals=None):
    """Markdown tables rendered straight from the result files.

    The effect table is built from the same registry values the prose uses,
    not recomputed from the CSV. Two code paths producing the same statistic
    can drift apart silently; audit_renderer.py fails on a registry key that
    nothing reads, which is what exposed the duplication.
    """
    v = vals if vals is not None else collect()
    t = {}
    rows = []
    for key, lab in (("all", "all one-to-one orthologs, no rhythm threshold"),
                     ("either", "either dataset q < 0.05"),
                     ("both", "both datasets q < 0.05"),
                     ("top200", "top 200 by combined rank")):
        if f"within_{key}" not in v:
            continue
        rows.append(f"| {lab} | {v[f'within_{key}']:.1f}% | "
                    f"{v[f'cross_{key}_median']:.1f}% | {v[f'cross_{key}_max']:.1f}% | "
                    f"{v[f'cross_{key}_sig']} of 10 |")
    t["effect_table"] = (
        "| ortholog selection | within strain | between lineages, median | max | "
        "pairs significant |\n|---|---|---|---|---|\n" + "\n".join(rows))

    # NOTE: a "power_table" generator built from n12_power_grid.csv used to sit
    # here. n12 is the 40-trial pilot; the template stopped referencing the
    # generator when the 400-trial surface replaced it, but the generator stayed,
    # still reading the pilot. Nothing failed, because the renderer audit checked
    # for unused registry KEYS and not for unused table GENERATORS. It now checks
    # both, and this generator is gone rather than left as a second power table
    # that no longer agrees with the one the manuscript prints.

    tn = _csv("n13_three_nulls_1000.csv")
    ci = "FPR_timepoint_ci95" in tn.columns
    cols = ("FPR_uniform", "FPR_timepoint", "FPR_marginal")
    rows = []
    for _, r in tn.iterrows():
        lab = ("none (uniform truth)" if r.frac_in_6h_window == 0
               else f"{r.frac_in_6h_window:.0%} in a 6 h window")
        mid = f"**{r[cols[1]]:.2f}**" + (f" ({r.FPR_timepoint_ci95})" if ci else "")
        rows.append(f"| {int(r.n_orthologs)} | {lab} | {r[cols[0]]:.2f} | {mid} | {r[cols[2]]:.2f} |")
    t["three_null_table"] = (
        "| n orthologs | phase concentration | uniform-phase | time-point permutation"
        + (" (95% CI)" if ci else "") + " | ortholog-identity |\n|---|---|---|---|---|\n"
        + "\n".join(rows))
    rb = _csv("n3_rbh_orthologs.csv")
    ko = _csv("n2_threshold_sensitivity.csv")
    ko = ko[ko.rule == "both q<0.05"]
    lab = {"SE7942_ito_WT": "*S. elongatus*", "SE7942_vijayan": "*S. elongatus* (study 2)",
           "SYN6803": "*Synechocystis*", "TERY": "*Trichodesmium*",
           "MED4": "*Prochlorococcus*", "NPUN": "*Nostoc*"}
    rows = []
    for _, r in rb.iterrows():
        within = {r.a, r.b} == WITHIN
        kv = ko[((ko.a == r.a) & (ko.b == r.b)) | ((ko.a == r.b) & (ko.b == r.a))]
        kos = f"{kv.effect_pct.iloc[0]:.1f}%" if len(kv) else "n/a"
        name = ("within strain, 7942 ×2" if within else
                f"{lab.get(r.a, r.a)} × {lab.get(r.b, r.b)}")
        rows.append(f"| {name} | {int(r.n_rbh)} | "
                    f"{r['effect_pct | all orthologs']:.1f}% | "
                    f"{r['effect_pct | both q<0.05']:.1f}% (p={r['p | both q<0.05']:.3f}) | {kos} |")
    ep = _csv("n11_effect_vs_pileup.csv")
    ep = ep[(ep.rule == "both q<0.05") &
            (~ep.apply(lambda r: {r.a, r.b} == WITHIN, axis=1))].sort_values(
                "effect_pct", ascending=False)
    rows = [f"| {lab.get(r.a, r.a)} × {lab.get(r.b, r.b)} | {int(r.n)} | "
            f"{r.effect_pct:+.1f} | {r.pile_mean:.2f} |" for _, r in ep.iterrows()]
    t["pair_table"] = ("| pair | orthologs | effect % | mean phase concentration |\n"
                       "|---|---|---|---|\n" + "\n".join(rows))

    # --- functional stratification, both BRITE levels
    lab2 = {"SE7942_ito_WT": "*S. elongatus*", "SE7942_vijayan": "*S. elongatus* (study 2)",
            "SYN6803": "*Synechocystis*", "TERY": "*Trichodesmium*",
            "MED4": "*Prochlorococcus*", "NPUN": "*Nostoc*"}
    for lvl in ("A", "B"):
        f = OUT / f"n18_stratified_closure_{lvl}.csv"
        if not f.exists():
            continue
        s = pd.read_csv(f)
        rows = []
        for _, r in s.iterrows():
            nm = ("**within strain, 7942 ×2**" if {r.a, r.b} == WITHIN
                  else f"{lab2.get(r.a, r.a)} × {lab2.get(r.b, r.b)}")
            rows.append(f"| {nm} | {int(r.n)} | {int(r.n_classes)} | {r.median_size:.0f} | "
                        f"{r.all_free_effect:.1f}% | {r.all_strat_effect:.1f}% (p={r.all_strat_p:.3f}) | "
                        f"{int(r.big_classes_only_n)} | {r.big_classes_only_strat_effect:.1f}% "
                        f"(p={r.big_classes_only_strat_p:.3f}) |")
        t[f"strat_table_{lvl}"] = (
            "| pair | orthologs | classes | median class size | free null | stratified |"
            " n in classes ≥10 | stratified, classes ≥10 |\n|---|---|---|---|---|---|---|---|\n"
            + "\n".join(rows))

    # --- design-matched null calibration
    f = OUT / "n16_design_matched_nulls.csv"
    if f.exists():
        d = pd.read_csv(f)
        for tag, sub in (("ours", d[d.null == "marginal"]), ("published", d[d.null == "timepoint"])):
            rows = [f"| {r.design} | {int(r.n_timepoints)} | {int(r.n_samples)} | "
                    f"{r.frac_in_6h_window:.0%} | {r.FPR:.3f} | {r.FPR_ci95} |"
                    for _, r in sub.iterrows()]
            t[f"design_table_{tag}"] = (
                "| design | timepoints | samples | phase concentration | FPR | 95% CI |\n"
                "|---|---|---|---|---|---|\n" + "\n".join(rows))

    # --- power surface (merged grid with Monte Carlo intervals)
    parts = [OUT / "n14_power_surface.csv", OUT / "n17_power_small_fractions.csv"]
    got = [pd.read_csv(x) for x in parts if x.exists()]
    if got:
        ps = pd.concat(got, ignore_index=True).sort_values(
            ["true_conserved_frac", "jitter_sd_h", "n_orthologs"])
        rows = []
        for c in sorted(ps.true_conserved_frac.unique()):
            for s in sorted(ps.jitter_sd_h.unique()):
                cells = []
                for n in (150, 500, 1500):
                    q = ps[(ps.true_conserved_frac == c) & (ps.jitter_sd_h == s)
                           & (ps.n_orthologs == n)]
                    cells.append(f"{q.power.iloc[0]:.2f} ({q.power_ci95.iloc[0]})"
                                 if len(q) else "—")
                rows.append(f"| {c:.0%} | {s:.0f} h | " + " | ".join(cells) + " |")
        t["power_surface_table"] = (
            "| true conserved fraction | jitter SD | n=150 | n=500 | n=1500 |\n"
            "|---|---|---|---|---|\n" + "\n".join(rows))

    # ---- design-matched calibration, reported as two separate tables
    #
    # (A) is the calibration this manuscript's own statistics rest on: the
    # primary null under the five sampling designs actually analysed.
    # (B) is a design-matched emulation of three published sampling designs
    # under a time-point-permutation null. It shares those designs' sampling
    # geometry only - not their organisms, their data, or their rhythm caller -
    # so it says what that null's size is under those designs when the truth is
    # zero by construction, and nothing about any published result.
    dm = _csv("n16_design_matched_nulls.csv")
    if dm is not None and len(dm):
        for key, which, caption in (("design_table_ours", "marginal", None),
                                    ("design_table_published", "timepoint", None)):
            d = dm[dm.null == which].sort_values(["design", "frac_in_6h_window"])
            if not len(d):
                continue
            fracs = sorted(d.frac_in_6h_window.unique())
            head = ("| sampling design | time points | samples | "
                    + " | ".join(("uniform truth" if f == 0 else f"{f:.0%} in 6 h")
                                 for f in fracs) + " |")
            rule = "|---" * (3 + len(fracs)) + "|"
            rows = []
            for name in d.design.unique():
                g = d[d.design == name]
                cells = []
                for f in fracs:
                    q = g[g.frac_in_6h_window == f]
                    cells.append(f"{q.FPR.iloc[0]:.3f} ({q.FPR_ci95.iloc[0]})"
                                 if len(q) else "—")
                r0 = g.iloc[0]
                rows.append(f"| {name} | {int(r0.n_timepoints)} | "
                            f"{int(r0.n_samples)} | " + " | ".join(cells) + " |")
            t[key] = "\n".join([head, rule] + rows)

    t["rbh_table"] = ("| pair | RBH orthologs | effect, all orthologs | effect, both q<0.05 | "
                      "same pair by KO |\n|---|---|---|---|---|\n" + "\n".join(rows))
    return t


if __name__ == "__main__":
    v = collect()
    (OUT / "numbers.json").write_text(json.dumps(v, indent=2, sort_keys=True))
    print(json.dumps(v, indent=2, sort_keys=True))
