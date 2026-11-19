"""Functional-stratified ortholog permutation, as a sensitivity analysis.

The primary null permutes ortholog identity freely, which assumes orthologs are
exchangeable once each dataset's marginal phase distribution is preserved. They
are not functionally interchangeable: if two lineages share timing among, say,
photosynthesis genes but not among division genes, free permutation destroys
that structure.

This null permutes identity only WITHIN a broad functional class, taken from the
KEGG BRITE ko00001 hierarchy - an annotation built independently of any phase
measurement in these datasets. It answers a more conditional question than the
primary null:

    is ortholog phase correspondence greater than expected given BOTH the
    marginal phase distributions AND broad functional composition?

A weaker result under this null does not invalidate the primary result. It would
indicate that part of the apparent correspondence is carried by functional
composition rather than by individual ortholog identity, which is a different
and also interesting statement. Both are reported.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import orthology as O
from analyze import OUT
from paths import KEGG as CACHE
from null_comparison import prepare, phases, QCUT
from three_nulls_fast import stat_many, PERIOD

MIN_CLASS = 10          # classes smaller than this are pooled, see brite_classes()


def brite_classes(level="B", multi="first"):
    """KO -> broad functional class from the KEGG BRITE ko00001 hierarchy.

    Level A is ~6 top categories (Metabolism, Genetic Information Processing, ...);
    level B is ~45 subcategories (Carbohydrate metabolism, Translation, ...).
    Neither uses expression data of any kind.

    **Multiply-assigned KOs.** BRITE is not a partition: 33% of KOs sit under
    more than one level-A category and 39% under more than one level-B
    category, up to 26 classes for one KO. Stratification needs one class per
    ortholog, so a rule is required, and the rule is arbitrary. `multi`
    selects it - "first" takes the first occurrence in the file (the default,
    and what the reported analysis uses), "last" the final occurrence, "alpha"
    the alphabetically first class name. closure_report compares the rules so
    a reader can see whether any conclusion depends on the choice.
    """
    txt = O._cached("https://rest.kegg.jp/get/br:ko00001", CACHE / "brite_ko00001.keg")
    out, a, b = {}, None, None
    for line in txt.splitlines():
        if not line or line[0] not in "ABCD":
            continue
        tag, rest = line[0], line[1:].strip()
        rest = rest.replace("<b>", "").replace("</b>", "").strip()
        if tag == "A":
            a = rest
        elif tag == "B":
            b = rest
        elif tag == "D":
            ko = rest.split()[0]
            if ko.startswith("K"):
                cls = a if level == "A" else b
                if multi == "first":
                    out.setdefault(ko, cls)
                elif multi == "last":
                    out[ko] = cls
                elif multi == "alpha":
                    prev = out.get(ko)
                    out[ko] = cls if prev is None or cls < prev else prev
                else:
                    raise ValueError(f"unknown multi rule: {multi}")
    return out


def stratified_perm(y, classes, rng):
    """permute y within each class only."""
    out = y.copy()
    for cls in np.unique(classes):
        idx = np.flatnonzero(classes == cls)
        if len(idx) > 1:
            out[idx] = rng.permutation(y[idx])
    return out


def compare(a, b, level="B", n_perm=1000, seed=0):
    ea, ma = prepare(a)
    eb, mb = prepare(b)
    shared = ea.index.intersection(eb.index)
    if len(shared) < 50:
        return None
    pa, qa = phases(ea.loc[shared], ma)
    pb, qb = phases(eb.loc[shared], mb)
    keep = ((qa < QCUT) & (qb < QCUT))
    keep = keep[keep].index
    if len(keep) < 30:
        return None
    kmap = brite_classes(level)
    cls = np.array([kmap.get(k, "unclassified") for k in keep])
    # pool classes below MIN_CLASS so the permutation is not degenerate
    counts = pd.Series(cls).value_counts()
    small = set(counts[counts < MIN_CLASS].index)
    cls_pooled = np.array(["pooled-small" if c in small else c for c in cls])
    x, y = pa[keep].values, pb[keep].values
    obs = stat_many(x, y[None, :])[0]
    rng = np.random.default_rng(seed)
    free = stat_many(x, np.array([rng.permutation(y) for _ in range(n_perm)]))
    strat = stat_many(x, np.array([stratified_perm(y, cls_pooled, rng)
                                   for _ in range(n_perm)]))
    sizes = pd.Series(cls_pooled).value_counts()
    permutable = int(sizes[sizes > 1].sum())

    def pack(nl, tag):
        nm = float(nl.mean())
        return {f"effect_{tag}": round((nm - obs) / nm * 100, 1),
                f"p_{tag}": round(float((nl <= obs).mean()), 4),
                f"nullmean_{tag}": round(nm, 3)}
    return dict(a=a, b=b, n_orthologs=len(keep), brite_level=level,
                n_classes=int(sizes.size), n_permutable=permutable,
                n_pooled_small=int(sizes.get("pooled-small", 0)),
                n_unclassified=int((cls == "unclassified").sum()),
                largest_class=int(sizes.max()), obs_mean_dphase=round(obs, 3),
                **pack(free, "free"), **pack(strat, "stratified"))


def main():
    pairs = [("SE7942_ito_WT", "SE7942_vijayan"), ("SE7942_ito_WT", "SYN6803"),
             ("SE7942_ito_WT", "TERY"), ("SE7942_ito_WT", "MED4"),
             ("SE7942_ito_WT", "NPUN"), ("TERY", "NPUN"), ("SYN6803", "NPUN")]
    rows = []
    for lvl in ("A", "B"):
        for a, b in pairs:
            r = compare(a, b, lvl)
            if r:
                rows.append(r)
                print({k: r[k] for k in ("a", "b", "brite_level", "n_orthologs",
                                         "n_classes", "n_permutable", "effect_free",
                                         "p_free", "effect_stratified", "p_stratified")},
                      flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "n15_stratified_null.csv", index=False)
    print("\n== free vs functionally stratified ortholog permutation")
    print(t[["a", "b", "brite_level", "n_orthologs", "n_classes", "n_permutable",
             "effect_free", "p_free", "effect_stratified", "p_stratified"]]
          .to_string(index=False))


def closure_report(level="B", n_perm=1000, seed=0, multi="first"):
    """Stratum sizes, small-stratum sensitivity, and the exchangeability note.

    The stratified null replaces one exchangeability assumption with another:
    orthologs are taken to be exchangeable WITHIN a broad functional class. That
    is weaker than full exchangeability but it is still an assumption, and broad
    BRITE classes can hide real differences in expression level and rhythm
    amplitude. These numbers let a reader judge it.
    """
    pairs = [("SE7942_ito_WT", "SE7942_vijayan"), ("SE7942_ito_WT", "SYN6803"),
             ("SE7942_ito_WT", "TERY"), ("SE7942_ito_WT", "MED4"),
             ("SE7942_ito_WT", "NPUN"), ("TERY", "NPUN"), ("SYN6803", "NPUN")]
    kmap = brite_classes(level, multi=multi)
    rows = []
    for a, b in pairs:
        ea, ma = prepare(a); eb, mb = prepare(b)
        shared = ea.index.intersection(eb.index)
        pa, qa = phases(ea.loc[shared], ma)
        pb, qb = phases(eb.loc[shared], mb)
        keep = ((qa < QCUT) & (qb < QCUT)); keep = keep[keep].index
        if len(keep) < 30:
            continue
        cls = pd.Series([kmap.get(k, "unclassified") for k in keep])
        sizes = cls.value_counts()
        x, y = pa[keep].values, pb[keep].values
        rng = np.random.default_rng(seed)
        obs = stat_many(x, y[None, :])[0]
        out = dict(a=a, b=b, level=level, multi_rule=multi, n_permutations=n_perm,
                   n=len(keep), n_classes=int(sizes.size),
                   largest=int(sizes.max()), median_size=float(sizes.median()),
                   n_singleton_classes=int((sizes == 1).sum()),
                   n_in_classes_lt10=int(sizes[sizes < 10].sum()),
                   n_unclassified=int(sizes.get("unclassified", 0)))
        # sensitivity: drop orthologs in classes smaller than 10 entirely
        big = cls.isin(sizes[sizes >= 10].index).values
        for tag, mask in (("all", np.ones(len(keep), bool)), ("big_classes_only", big)):
            if mask.sum() < 30:
                continue
            xs, ys = x[mask], y[mask]
            cs = cls.values[mask]
            o = stat_many(xs, ys[None, :])[0]
            free = stat_many(xs, np.array([rng.permutation(ys) for _ in range(n_perm)]))
            strat = stat_many(xs, np.array([stratified_perm(ys, cs, rng)
                                            for _ in range(n_perm)]))
            for nm, nl in (("free", free), ("strat", strat)):
                m = float(nl.mean())
                out[f"{tag}_{nm}_effect"] = round((m - o) / m * 100, 1)
                out[f"{tag}_{nm}_p"] = round(float((nl <= o).mean()), 4)
            out[f"{tag}_n"] = int(mask.sum())
        rows.append(out)
        print(out, flush=True)
    t = pd.DataFrame(rows)
    suffix = level if multi == "first" else f"{level}_{multi}"
    t.to_csv(OUT / f"n18_stratified_closure_{suffix}.csv", index=False)
    return t


def multiassign_report():
    """How far BRITE departs from being a partition, per level.

    Stratification requires one class per ortholog, but a third of KOs carry
    several. The numbers below are what makes the choice of rule a reportable
    analytical decision rather than an implementation detail.
    """
    from collections import defaultdict
    txt = O._cached("https://rest.kegg.jp/get/br:ko00001", CACHE / "brite_ko00001.keg")
    rows = []
    for level in ("A", "B"):
        seen = defaultdict(set)
        a = b = None
        for line in txt.splitlines():
            if not line or line[0] not in "ABCD":
                continue
            tag, rest = line[0], line[1:].strip()
            rest = rest.replace("<b>", "").replace("</b>", "").strip()
            if tag == "A":
                a = rest
            elif tag == "B":
                b = rest
            elif tag == "D":
                ko = rest.split()[0]
                if ko.startswith("K"):
                    seen[ko].add(a if level == "A" else b)
        multi = sum(1 for v in seen.values() if len(v) > 1)
        rows.append(dict(brite_level=level, n_kos=len(seen),
                         n_classes=len({c for v in seen.values() for c in v}),
                         n_multiassigned=multi,
                         pct_multiassigned=round(100 * multi / len(seen), 1),
                         max_classes_per_ko=max(len(v) for v in seen.values())))
        print(rows[-1], flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "n19_brite_multiassign.csv", index=False)
    return t


if __name__ == "__main__":
    main()
    for lvl in ("A", "B"):
        print(f"\n== stratified-null closure report, BRITE level {lvl}")
        closure_report(lvl)
    print("\n== how far BRITE departs from a partition")
    multiassign_report()
    # sensitivity to the arbitrary rule for multiply-assigned KOs
    for lvl in ("A", "B"):
        print(f"\n== level {lvl} again, multiply-assigned KOs by last occurrence")
        closure_report(lvl, multi="last")
