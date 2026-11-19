"""Repeat the key comparison on reciprocal-best-hit orthologs instead of KEGG KOs.

If the within-strain vs cross-lineage gap is an artifact of how orthology was
assigned, swapping KEGG Orthology for reciprocal best hits computed directly
between the two proteomes should change it.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import load, rhythm, orthology as O
from analyze import OUT, phase_diff
from null_comparison import GEO, best_shift, QCUT

SHIFTS = np.arange(-12, 12.5, 1.0)


def best_shift_many(x, Y, period=24.0, block=64):
    """Vectorised best-shift mean |dphase| for many permuted vectors at once."""
    out = np.empty(len(Y))
    for i in range(0, len(Y), block):
        Yb = Y[i:i + block]
        d = (x[None, None, :] - Yb[:, None, :] - SHIFTS[None, :, None]) % period
        d = np.minimum(d, period - d)
        out[i:i + block] = d.mean(axis=2).min(axis=1)
    return out
import rbh as RBH

GENOME_OF = {"SE7942_ito_WT": "SE7942", "SE7942_vijayan": "SE7942",
             "SYN6803": "SYN6803", "TERY": "TERY", "MED4": "MED4", "NPUN": "NPUN"}
BRIDGE = {"MED4": "MED4.gff.gz"}


def expr_by_original_tag(did):
    """Expression indexed by the ORIGINAL locus tag, which is what the RBH maps use."""
    expr, meta = load.LOADERS[did](GEO)
    if did in BRIDGE:
        b = O.locus_tag_bridge(GEO / BRIDGE[did])
        idx = [b.get(g, "") for g in expr.index]
        expr = expr.loc[[i != "" for i in idx]]
        expr.index = [i for i in idx if i != ""]
    expr = expr[~expr.index.duplicated()]
    return expr, meta


def compare(a, b, n_perm=1000, seed=0):
    ga, gb = GENOME_OF[a], GENOME_OF[b]
    if ga == gb:
        pairs = None                       # same genome: orthology is identity
    else:
        try:
            pairs = RBH.rbh(ga, gb, cache=True)
        except Exception:
            return None
        if not pairs:
            return None
    ea, ma = expr_by_original_tag(a)
    eb, mb = expr_by_original_tag(b)
    if pairs is None:
        shared = ea.index.intersection(eb.index)
        xa, xb = list(shared), list(shared)
    else:
        xa = [k for k, v in pairs.items() if k in ea.index and v in eb.index]
        xb = [pairs[k] for k in xa]
    if len(xa) < 30:
        return None
    ra = rhythm.cosinor(ea.loc[xa], ma)
    rb = rhythm.cosinor(eb.loc[xb], mb)
    pa, pb = ra["phase_h"].values, rb["phase_h"].values
    qa, qb = ra["q"].values, rb["q"].values
    rng = np.random.default_rng(seed)
    out = {}
    for name, mask in (("all orthologs", np.ones(len(pa), bool)),
                       ("both q<0.05", (qa < QCUT) & (qb < QCUT))):
        if mask.sum() < 20:
            continue
        x, y = pa[mask], pb[mask]
        obs = float(best_shift_many(x, y[None, :])[0])
        null = best_shift_many(x, np.array([rng.permutation(y) for _ in range(n_perm)]))
        nm = float(null.mean())
        out[name] = dict(n=int(mask.sum()), obs=round(obs, 2), null=round(nm, 2),
                         effect_pct=round((nm - obs) / nm * 100, 1),
                         p=round(float((null <= obs).mean()), 4))
    return dict(a=a, b=b, n_rbh=(len(pairs) if pairs else len(xa)), **{
        f"{k} | {m}": v for m, d in out.items() for k, v in d.items()})


def main():
    pairs = [("SE7942_ito_WT", "SE7942_vijayan"), ("SE7942_ito_WT", "SYN6803"),
             ("SE7942_ito_WT", "TERY"), ("SE7942_ito_WT", "MED4"), ("TERY", "MED4")]
    rows = [r for r in (compare(a, b) for a, b in pairs) if r]
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "n3_rbh_orthologs.csv", index=False)
    print(t.to_string(index=False))


def write_params_record():
    """Persist the RBH parameters and the checksum of every input proteome."""
    import hashlib
    import json
    from rbh import RBH_PARAMS
    from paths import PROTEOMES
    rec = {"params": RBH_PARAMS, "proteomes": {}}
    for f in sorted(PROTEOMES.glob("*")):
        if f.is_file():
            h = hashlib.sha256(f.read_bytes()).hexdigest()
            rec["proteomes"][f.name] = {"sha256": h, "bytes": f.stat().st_size}
    out = OUT / "n3_rbh_params.json"
    out.write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n")
    print(f"wrote {out.name}: {len(rec['proteomes'])} proteomes recorded")
    return rec


if __name__ == "__main__":
    main()
    write_params_record()
