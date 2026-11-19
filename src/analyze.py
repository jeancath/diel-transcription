"""Analyses on top of the per-gene rhythm calls.

Q1  Is the circadian transcriptional program reproducible *within* the model
    strain, across independent studies?  (the within-species ceiling)
Q2  Is it conserved *across* lineages at the level of orthologous genes?
Q3  Are the phases of the conserved-rhythmic orthologs conserved?
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import orthology as O

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
from paths import KEGG as CACHE
QCUT = 0.05
CLOCK_KO = {"K08480": "kaiA", "K08481": "kaiB", "K08482": "kaiC",
            "K08479": "sasA", "K10697": "rpaA", "K29143": "cirA"}


def circ_mean(h, period=24.0):
    a = np.asarray(h, float) * 2 * np.pi / period
    return (np.arctan2(np.sin(a).mean(), np.cos(a).mean()) * period / (2 * np.pi)) % period


def circ_sd(h, period=24.0):
    """Circular SD in hours; 0 = identical phase, ~6.9 h = uniform on a 24 h circle."""
    a = np.asarray(h, float) * 2 * np.pi / period
    R = np.hypot(np.sin(a).mean(), np.cos(a).mean())
    R = min(max(R, 1e-12), 1 - 1e-12)
    return float(np.sqrt(-2 * np.log(R)) * period / (2 * np.pi))


def phase_diff(a, b, period=24.0):
    d = (np.asarray(a, float) - np.asarray(b, float)) % period
    return np.where(d > period / 2, d - period, d)


def load():
    df = pd.read_csv(OUT / "rhythms_per_gene.csv.gz")
    df["rhythmic"] = df["q"] < QCUT
    return df


# ------------------------------------------------------------------ Q1
def within_species(df):
    a = df[df.dataset == "SE7942_ito_WT"].set_index("gene")
    b = df[df.dataset == "SE7942_vijayan"].set_index("gene")
    common = a.index.intersection(b.index)
    a, b = a.loc[common], b.loc[common]
    ra, rb = set(a.index[a.rhythmic]), set(b.index[b.rhythmic])
    both = ra & rb
    jac = len(both) / len(ra | rb)
    # phase agreement, restricted to genes both studies call rhythmic
    d = phase_diff(a.loc[sorted(both), "phase_h"], b.loc[sorted(both), "phase_h"])
    within2 = float((np.abs(d) <= 2).mean())
    within4 = float((np.abs(d) <= 4).mean())
    # expectation if the two studies agreed only by chance
    exp_overlap = len(ra) * len(rb) / len(common)
    rows = [dict(metric="genes compared (both studies)", value=len(common)),
            dict(metric="rhythmic, Ito 2009 (LL, GPL8012)", value=len(ra)),
            dict(metric="rhythmic, Vijayan 2009 (LL, GPL9534)", value=len(rb)),
            dict(metric="rhythmic in both", value=len(both)),
            dict(metric="expected overlap if independent", value=round(exp_overlap, 1)),
            dict(metric="Jaccard of rhythmic sets", value=round(jac, 3)),
            dict(metric="replication rate of Ito calls in Vijayan",
                 value=round(len(both) / len(ra), 3)),
            dict(metric="phase agreement within 2 h", value=round(within2, 3)),
            dict(metric="phase agreement within 4 h", value=round(within4, 3)),
            dict(metric="median |phase difference| (h)",
                 value=round(float(np.median(np.abs(d))), 2)),
            dict(metric="circular SD of phase difference (h)", value=round(circ_sd(d), 2))]
    return pd.DataFrame(rows), pd.DataFrame({"gene": sorted(both), "phase_diff_h": d})


# ------------------------------------------------------------------ Q2 / Q3
def cross_lineage(df, datasets):
    sub = df[df.dataset.isin(datasets) & df.ko.notna() & (df.ko != "")].copy()
    # one row per (dataset, KO): keep the most significant gene of a KO family
    sub = sub.sort_values("q").drop_duplicates(["dataset", "ko"])
    piv_r = sub.pivot(index="ko", columns="dataset", values="rhythmic")
    piv_p = sub.pivot(index="ko", columns="dataset", values="phase_h")
    piv_a = sub.pivot(index="ko", columns="dataset", values="amp_log2")
    shared = piv_r.dropna()
    n_ds = shared.shape[1]
    tab = pd.DataFrame(index=shared.index)
    tab["n_rhythmic"] = shared.sum(axis=1).astype(int)
    tab["phase_circ_sd_h"] = [circ_sd(piv_p.loc[k, shared.columns].values) for k in shared.index]
    tab["phase_circ_mean_h"] = [circ_mean(piv_p.loc[k, shared.columns].values) for k in shared.index]
    # phase SD computed only where the KO is actually called rhythmic
    sds = []
    for k in shared.index:
        ph = piv_p.loc[k, shared.columns][shared.loc[k].astype(bool).values]
        sds.append(circ_sd(ph.values) if len(ph) >= 3 else np.nan)
        tab.loc[k, "mean_amp_log2"] = float(np.nanmean(piv_a.loc[k, shared.columns].values))
    tab["phase_sd_rhythmic_only_h"] = sds
    names = O.ko_names(CACHE)
    tab["ko_name"] = [names.get(k, "")[:70] for k in tab.index]
    return tab, shared, piv_p, n_ds


def binomial_expectation(shared):
    """Per-dataset rhythmic rates -> expected count of KOs rhythmic in k datasets
    if rhythmicity were independent across lineages."""
    rates = shared.mean(axis=0).values
    n = shared.shape[0]
    # Poisson-binomial by convolution
    dist = np.array([1.0])
    for r in rates:
        dist = np.convolve(dist, [1 - r, r])
    return pd.DataFrame({"k_datasets": np.arange(len(dist)),
                         "expected_KOs": (dist * n).round(1),
                         "observed_KOs": [int((shared.sum(axis=1) == k).sum())
                                          for k in range(len(dist))]})


def main():
    df = load()
    OUT.mkdir(exist_ok=True)

    # ---- positive control: the clock genes themselves
    pc = df[df.ko.isin(CLOCK_KO)][["dataset", "lineage", "light", "gene", "ko",
                                   "amp_log2", "phase_h", "q"]].copy()
    pc["clock_gene"] = pc.ko.map(CLOCK_KO)
    pc.sort_values(["clock_gene", "dataset"]).to_csv(OUT / "clock_gene_control.csv", index=False)

    # ---- Q1
    w, pd_genes = within_species(df)
    w.to_csv(OUT / "q1_within_species_reproducibility.csv", index=False)
    print("== Q1  Synechococcus elongatus PCC 7942, two independent LL studies")
    print(w.to_string(index=False))

    # ---- Q2/Q3: real lineages only, drop controls and the underpowered LL sets
    test_sets = ["SE7942_ito_WT", "SYN6803", "TERY", "MED4", "NPUN"]
    tab, shared, piv_p, n_ds = cross_lineage(df, test_sets)
    tab.sort_values(["n_rhythmic", "phase_sd_rhythmic_only_h"],
                    ascending=[False, True]).to_csv(OUT / "q2_ko_conservation.csv")
    exp = binomial_expectation(shared)
    exp.to_csv(OUT / "q2_expected_vs_observed.csv", index=False)
    print(f"\n== Q2  {shared.shape[0]} KEGG orthologs measured in all {n_ds} lineages")
    print(exp.to_string(index=False))
    core = tab[tab.n_rhythmic == n_ds]
    print(f"\nrhythmic in all {n_ds}: {len(core)} KOs "
          f"({len(core)/shared.shape[0]:.1%} of shared orthologs)")
    print(f"median phase circular SD among those: "
          f"{core['phase_sd_rhythmic_only_h'].median():.2f} h "
          f"(uniform phase would be ~6.9 h)")
    print("\ntop conserved-rhythmic orthologs by phase agreement:")
    print(core.sort_values("phase_sd_rhythmic_only_h")
              .head(15)[["ko_name", "phase_circ_mean_h", "phase_sd_rhythmic_only_h"]]
              .to_string())
    # phase conservation as a function of how many lineages call it rhythmic
    g = (tab.dropna(subset=["phase_sd_rhythmic_only_h"])
            .groupby("n_rhythmic")["phase_sd_rhythmic_only_h"]
            .agg(["count", "median"]).round(2))
    print("\nphase circular SD by number of lineages calling the ortholog rhythmic:")
    print(g.to_string())
    g.to_csv(OUT / "q3_phase_sd_by_conservation.csv")


if __name__ == "__main__":
    main()
