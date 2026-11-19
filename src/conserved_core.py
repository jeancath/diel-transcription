"""What is in the weak conserved component?

The cross-lineage effect is small but non-zero. This asks which orthologs carry
it and what they do. Each lineage is first aligned to the reference strain by
the single global time shift that minimises mean |phase difference| over all its
one-to-one orthologs - removing the arbitrary per-study time origin - and the
dispersion of each ortholog's phase across lineages is then measured on the
aligned phases. Functional enrichment is hypergeometric against the background
of orthologs measured in every lineage, which is the right background: it asks
what distinguishes the phase-conserved orthologs from the comparable set, not
from the genome.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import orthology as O
from analyze import OUT, circ_sd, circ_mean, phase_diff
from null_comparison import prepare, phases, best_shift, CACHE, QCUT

REF = "SE7942_ito_WT"
LINEAGES = [REF, "SYN6803", "TERY", "MED4", "NPUN"]


def aligned_phase_table():
    ph, qs = {}, {}
    for d in LINEAGES:
        e, m = prepare(d)
        p, q = phases(e, m)
        ph[d], qs[d] = p, q
    ref = ph[REF]
    shifts = {REF: 0.0}
    for d in LINEAGES[1:]:
        common = ref.index.intersection(ph[d].index)
        _, s = best_shift(ref[common].values, ph[d][common].values)
        shifts[d] = s
    P = pd.DataFrame({d: (ph[d] + shifts[d]) % 24 for d in LINEAGES})
    Q = pd.DataFrame({d: qs[d] for d in LINEAGES})
    return P, Q, shifts


def main():
    P, Q, shifts = aligned_phase_table()
    print("global alignment shift applied to each lineage (h, relative to PCC 7942):")
    print("  " + ", ".join(f"{d}: {s:+.0f}" for d, s in shifts.items()))

    shared = P.dropna().index
    P, Q = P.loc[shared], Q.loc[shared]
    rhythmic = (Q < QCUT)
    tab = pd.DataFrame(index=shared)
    tab["n_rhythmic"] = rhythmic.sum(axis=1)
    tab["phase_sd_h"] = [circ_sd(P.loc[k][rhythmic.loc[k]].values)
                         if rhythmic.loc[k].sum() >= 3 else np.nan for k in shared]
    tab["phase_mean_h"] = [circ_mean(P.loc[k][rhythmic.loc[k]].values)
                           if rhythmic.loc[k].sum() >= 3 else np.nan for k in shared]
    names = O.ko_names(CACHE)
    tab["ko_name"] = [names.get(k, "")[:72] for k in shared]
    tab = tab.sort_values("phase_sd_h")
    tab.to_csv(OUT / "n5_conserved_core.csv")
    print(f"\n{len(shared)} one-to-one orthologs measured in all {len(LINEAGES)} lineages; "
          f"{int((tab.n_rhythmic >= 4).sum())} rhythmic in >=4")

    core = tab[(tab.n_rhythmic >= 4) & (tab.phase_sd_h < tab.phase_sd_h.quantile(.25))]
    print(f"\n== the {len(core)} most phase-conserved orthologs (aligned phases)")
    print(core[["n_rhythmic", "phase_mean_h", "phase_sd_h", "ko_name"]]
          .head(25).to_string())

    # --- pathway enrichment of the conserved core against the shared background
    txt = O._cached("https://rest.kegg.jp/link/pathway/ko", CACHE / "ko_pathway.tsv")
    pw = {}
    for line in txt.splitlines():
        if "\t" in line:
            k, p = line.split("\t")[:2]
            if p.startswith("path:map"):
                pw.setdefault(k.replace("ko:", ""), set()).add(p.replace("path:", ""))
    pnames = {}
    for line in O._cached("https://rest.kegg.jp/list/pathway", CACHE / "pathway_names.tsv").splitlines():
        if "\t" in line:
            a, b = line.split("\t", 1)
            pnames[a.replace("path:", "")] = b
    bg = set(shared)
    cs = set(core.index)
    rows = []
    allp = set().union(*[pw.get(k, set()) for k in bg]) if bg else set()
    for p in allp:
        inp = {k for k in bg if p in pw.get(k, set())}
        hit = inp & cs
        if len(hit) < 3:
            continue
        pval = stats.hypergeom.sf(len(hit) - 1, len(bg), len(inp), len(cs))
        rows.append(dict(pathway=p, name=pnames.get(p, "")[:58], n_bg=len(inp),
                         n_core=len(hit),
                         frac=round(len(hit) / len(inp), 2), p=pval))
    # stable, fully determined order: pandas' default sort is not stable, so
    # two pathways tied on p could swap places between identical builds
    en = pd.DataFrame(rows).sort_values(["p", "pathway"], kind="mergesort")
    if len(en):
        from rhythm import bh
        en["q"] = bh(en.p.values)
        en.to_csv(OUT / "n5_core_pathway_enrichment.csv", index=False)
        print("\n== pathways over-represented in the phase-conserved core")
        print(en.head(12).to_string(index=False, float_format="%.3g"))

    sd_null_report(P, Q)


def sd_null(P, Q, n_perm=500, seed=0, min_rhythmic=4, cutoff_from=3):
    """Is the phase-conserved core bigger than alignment alone would produce?

    The global alignment shifts are chosen to minimise mean |phase difference|,
    which makes the lineages' phase modes coincide; orthologs sitting near that
    shared mode will then have a low phase SD whether or not their timing is
    ortholog-specific. This null shuffles ortholog identity independently within
    each lineage - preserving every lineage's aligned phase distribution and the
    alignment itself - and recomputes the same per-ortholog SDs.

    The two criteria are separate on purpose, because main() defines the core
    with two different populations: membership requires a gene to be rhythmic
    in at least `min_rhythmic` lineages, while the dispersion cutoff is the
    lower quartile over every gene scorable at all (rhythmic in at least
    `cutoff_from`). Passing (4, 3) reproduces exactly the core the Results
    report; any other combination counts a different gene set and cannot be
    quoted as its control.

    Caveat that limits what a small p value here means: the alignment shifts
    were fitted using the true ortholog pairings and are held fixed across
    permutations, so coincidence induced by the fitted alignment is inside the
    observed statistic but not inside the null. The test is therefore not
    alignment-protected and is anti-conservative to an unquantified degree.
    """
    rng = np.random.default_rng(seed)
    rhythmic = (Q < QCUT)
    idx = P.index

    def score(Pv, Rv):
        sds = np.array([circ_sd(Pv[i][Rv[i]]) if Rv[i].sum() >= cutoff_from else np.nan
                        for i in range(len(idx))])
        enough = np.array([Rv[i].sum() >= min_rhythmic for i in range(len(idx))])
        return sds, enough

    Pv0, Rv0 = P.values, rhythmic.values
    sds, enough = score(Pv0, Rv0)
    cut = float(np.nanquantile(sds, .25))
    n_obs = int(np.nansum((sds < cut) & enough))
    counts = []
    for _ in range(n_perm):
        Pp = np.column_stack([rng.permutation(Pv0[:, j]) for j in range(Pv0.shape[1])])
        Rp = np.column_stack([rng.permutation(Rv0[:, j]) for j in range(Rv0.shape[1])])
        s2, e2 = score(Pp, Rp)
        counts.append(int(np.nansum((s2 < cut) & e2)))
    counts = np.array(counts)
    return dict(min_rhythmic=min_rhythmic, cutoff_from=cutoff_from,
                n_permutations=n_perm,
                sd_cutoff_h=round(cut, 2), observed_below_cutoff=n_obs,
                null_mean=round(float(counts.mean()), 1),
                null_p95=int(np.quantile(counts, .95)),
                p=round(float((counts >= n_obs).mean()), 4),
                enrichment=round(n_obs / counts.mean(), 2) if counts.mean() else np.nan)


def sd_null_report(P, Q, n_perm=500, seed=0):
    """Run the ortholog-identity control at both rhythmic-count criteria.

    The core reported in the Results requires a gene to be rhythmic in >= 4 of
    the 5 lineages, so the control that licenses a statement about that core
    must use the same criterion. An earlier version of this function used >= 3
    and was never wired into any make target, so the control reported in the
    manuscript counted a different gene set from the core it was controlling.
    Both criteria are reported here rather than silently choosing one.
    """
    rows = [sd_null(P, Q, n_perm=n_perm, seed=seed,
                    min_rhythmic=m, cutoff_from=c)
            for m, c in ((4, 3), (3, 3))]
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "n5_core_sd_null.csv", index=False)
    print("\n== ortholog-identity control on the phase-conserved core")
    print(t.to_string(index=False))
    return t


if __name__ == "__main__":
    main()
