"""Driver: call rhythms in every dataset, attach KEGG Orthology, write tables."""
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import load, rhythm, orthology as O
from datasets import DATASETS

from paths import GEO, KEGG as CACHE
OUT = Path(__file__).resolve().parents[1] / "results"
BRIDGE = {"MED4": "MED4.gff.gz", "HALO_05": "HAO.gff.gz", "HALO_25": "HAO.gff.gz"}
QCUT = 0.05


def main():
    CACHE.mkdir(exist_ok=True)
    OUT.mkdir(exist_ok=True)
    frames, summary = [], []
    for d in DATASETS:
        expr, meta = load.LOADERS[d.did](GEO)
        res = rhythm.cosinor(expr, meta, period=24.0)
        tau = rhythm.best_period(expr, meta)
        res_tau = rhythm.cosinor(expr, meta, period=tau)

        g2ko = O.gene_to_ko(d.kegg, CACHE)
        if d.did in BRIDGE:
            bridge = O.locus_tag_bridge(GEO / BRIDGE[d.did])
            res["kegg_gene"] = [bridge.get(g, "") for g in res.index]
        else:
            res["kegg_gene"] = res.index
        res["ko"] = [g2ko.get(g, "") for g in res["kegg_gene"]]
        res["dataset"] = d.did
        res["lineage"] = d.lineage
        res["light"] = d.light
        res["role"] = d.role
        res["q_tau"] = res_tau["q"].reindex(res.index).values
        res["phase_tau"] = res_tau["phase_h"].reindex(res.index).values
        res["tau_fit"] = tau
        res.index.name = "gene"
        frames.append(res.reset_index())

        n = len(res)
        rhy = (res["q"] < QCUT).sum()
        summary.append(dict(dataset=d.did, lineage=d.lineage, gse=d.gse, light=d.light,
                            platform=d.platform, role=d.role, n_samples=expr.shape[1],
                            n_timepoints=meta["time"].nunique(), n_genes=n,
                            n_rhythmic=int(rhy), frac_rhythmic=round(rhy / n, 4),
                            n_ko_annotated=int((res["ko"] != "").sum()),
                            tau_fit=tau,
                            n_rhythmic_tau=int((res["q_tau"] < QCUT).sum())))
        print(f"{d.did:20s} genes={n:5d} rhythmic(q<{QCUT})={rhy:5d} "
              f"({rhy/n:5.1%})  tau={tau:4.1f}h  rhythmic@tau={summary[-1]['n_rhythmic_tau']:5d}",
              flush=True)

    allres = pd.concat(frames, ignore_index=True)
    allres.to_csv(OUT / "rhythms_per_gene.csv.gz", index=False)
    pd.DataFrame(summary).to_csv(OUT / "dataset_summary.csv", index=False)
    print("\nwrote", OUT / "rhythms_per_gene.csv.gz")


if __name__ == "__main__":
    main()
