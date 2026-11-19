"""Robustness of the conservation result.

Two objections have to be answered before a chance-level co-rhythmicity result
means anything:

1. The per-dataset rhythmic sets are large (56-98% of measured genes), so a big
   raw overlap is uninformative. Fixed amplitude floors do not fix this either,
   because the microarray datasets compress amplitude relative to RNA-seq - a
   common absolute floor would silently discard the model strain. So rhythmicity
   is defined by RANK within each dataset (top X% by cosinor p), which equalises
   the marginal rate across platforms and makes the independence expectation
   exact.

2. Phases are not uniformly distributed in any dataset - transcripts pile up
   near dawn and dusk. Apparent cross-lineage phase agreement can therefore come
   entirely from similar marginal phase distributions rather than from anything
   ortholog-specific. The permutation null below shuffles phases among the same
   set of orthologs, preserving each dataset's marginal phase distribution, and
   so tests only the ortholog-specific component.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from analyze import OUT, circ_sd, phase_diff, load

TEST = ["SE7942_ito_WT", "SYN6803", "TERY", "MED4", "NPUN"]
WITHIN = ("SE7942_ito_WT", "SE7942_vijayan")


def rank_table(df, datasets, top_frac):
    """rhythmic = top `top_frac` of genes by cosinor p within each dataset."""
    sub = df[df.dataset.isin(datasets) & (df.ko != "") & df.ko.notna()].copy()
    sub = sub.sort_values("p").drop_duplicates(["dataset", "ko"])
    sub["rhythmic"] = False
    for d, g in sub.groupby("dataset"):
        n = max(1, int(round(top_frac * len(g))))
        sub.loc[g.nsmallest(n, "p").index, "rhythmic"] = True
    r = sub.pivot(index="ko", columns="dataset", values="rhythmic").dropna()
    p = sub.pivot(index="ko", columns="dataset", values="phase_h").loc[r.index]
    a = sub.pivot(index="ko", columns="dataset", values="amp_log2").loc[r.index]
    return r.astype(bool), p, a


def poisson_binomial(rates, n):
    dist = np.array([1.0])
    for x in rates:
        dist = np.convolve(dist, [1 - x, x])
    return dist * n


def rank_sweep(df, datasets=None, label=""):
    datasets = datasets or TEST
    rows = []
    for tf in (0.10, 0.25, 0.50):
        r, _, _ = rank_table(df, datasets, tf)
        k = r.shape[1]
        exp = poisson_binomial(r.mean(axis=0).values, r.shape[0])
        obs_all = int((r.sum(axis=1) == k).sum())
        rows.append(dict(set=label or "all", top_frac=tf, lineages=k,
                         shared_KOs=int(r.shape[0]),
                         rhythmic_in_all=obs_all,
                         expected_if_independent=round(float(exp[k]), 2),
                         obs_over_exp=(round(obs_all / exp[k], 2) if exp[k] > 0 else np.nan),
                         rhythmic_in_ge_k_minus_1=int((r.sum(axis=1) >= k - 1).sum()),
                         expected_ge_k_minus_1=round(float(exp[k - 1:].sum()), 2)))
    return pd.DataFrame(rows)


def phase_test(df, datasets, top_frac=0.25, n_perm=4000, seed=0, tol=3.0):
    """Ortholog-specific phase agreement against a marginal-preserving null.

    Each pair is ranked and intersected on its OWN shared orthologs, not on the
    orthologs shared by every lineage in the panel: restricting to the global
    intersection would throw away most of the comparable genes for no reason.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for i, a in enumerate(datasets):
        for b in datasets[i + 1:]:
            r, p, _ = rank_table(df, [a, b], top_frac)
            both = r[a] & r[b]
            if both.sum() < 15:
                continue
            x, y = p.loc[both, a].values, p.loc[both, b].values
            obs = float((np.abs(phase_diff(x, y)) <= tol).mean())
            null = np.array([(np.abs(phase_diff(x, rng.permutation(y))) <= tol).mean()
                             for _ in range(n_perm)])
            rows.append(dict(a=a, b=b, n_orthologs=int(both.sum()),
                             median_abs_dphase_h=round(float(np.median(np.abs(x - y) % 24)), 2),
                             obs_frac_within_3h=round(obs, 3),
                             null_mean=round(float(null.mean()), 3),
                             null_p95=round(float(np.quantile(null, 0.95)), 3),
                             enrichment=round(obs / null.mean(), 2) if null.mean() else np.nan,
                             p_perm=float((null >= obs).mean())))
    return pd.DataFrame(rows)


def phase_marginals(df):
    """Where do rhythmic transcripts peak in each dataset? (the dawn/dusk pile-up)"""
    rows = []
    for d, g in df[df.q < 0.05].groupby("dataset"):
        h = np.histogram(g.phase_h, bins=np.arange(0, 25, 3))[0]
        rows.append(dict(dataset=d, n=len(g),
                         **{f"ZT{b}-{b+3}": int(c) for b, c in zip(range(0, 24, 3), h)},
                         frac_in_busiest_6h=round(float(
                             max(h[i] + h[i + 1] for i in range(len(h) - 1)) / len(g)), 3)))
    return pd.DataFrame(rows)


def main():
    df = load()
    sweeps = [rank_sweep(df, TEST, "5 lineages"),
              rank_sweep(df, [d for d in TEST if d != "MED4"], "no Prochlorococcus"),
              rank_sweep(df, list(WITHIN), "within-strain (7942 x2)")]
    sw = pd.concat(sweeps, ignore_index=True)
    sw.to_csv(OUT / "r1_rank_sweep.csv", index=False)
    print("== Co-rhythmicity vs independence, rank-based rhythm definition")
    print(sw.to_string(index=False))

    pt = phase_test(df, TEST + ["SE7942_vijayan"])
    ld_only = ["SYN6803", "TERY", "MED4", "NPUN"]
    pt["comparison"] = np.where(pt[["a", "b"]].apply(lambda s: set(s) == set(WITHIN), axis=1),
                                "WITHIN-STRAIN", "cross-lineage")
    pt.to_csv(OUT / "r2_phase_test.csv", index=False)
    print("\n== Ortholog-specific phase agreement (top 25% by p, tolerance 3 h)")
    print(pt.sort_values("enrichment", ascending=False).to_string(index=False))

    sub = pt[pt.a.isin(ld_only) & pt.b.isin(ld_only)]
    print(f"\nLD-driven lineages only (same external cue, no LL/LD mismatch): "
          f"median enrichment {sub.enrichment.median():.2f}, "
          f"n pairs {len(sub)}, pairs with p<0.05: {(sub.p_perm < 0.05).sum()}")
    pm = phase_marginals(df)
    pm.to_csv(OUT / "r3_phase_marginals.csv", index=False)
    print("\n== Marginal phase distribution of rhythmic transcripts")
    print(pm.to_string(index=False))


# Three functions were removed from here: genome_similarity, decay and
# phase_test_offset_free. All three sat below the __main__ block, so no make
# target could ever reach them; none wrote a result file and no manuscript
# claim depended on them. They are in git history if the divergence-decay or
# offset-optimised pairwise analysis is ever wanted as a real analysis.


if __name__ == "__main__":
    main()
