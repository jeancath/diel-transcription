"""Loaders: GEO series matrices and supplementary count tables -> (log2 expression, sample metadata).

Sampling time is recovered from the GEO sample titles / characteristics for every
dataset. Where a dataset holds several experimental arms (clock mutants, drug
treatments, extended dark) only the arm named by the dataset id is kept.
"""
import gzip, io, re
import numpy as np
import pandas as pd


# ---------------------------------------------------------------- series matrix
def _series_matrix(path):
    """Return (table DataFrame, meta dict of header lines)."""
    meta, rows, in_table = {}, [], False
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                rows.append(line.rstrip("\n"))
            elif line.startswith("!"):
                k, _, v = line.rstrip("\n").partition("\t")
                meta.setdefault(k[1:], []).append([x.strip('"') for x in v.split("\t")])
    tab = pd.read_csv(io.StringIO("\n".join(rows)), sep="\t", index_col=0)
    tab.index = [str(i).strip('"') for i in tab.index]
    return tab, meta


def _titles(meta):
    return meta["Sample_title"][0]


def _probe_map(gpl_path, probe_col=0, orf_col=None, header_pick=None, pattern=None):
    """probe id -> gene/ORF id, from a GEO platform annotation file.

    `pattern` guards against column drift in these hand-made platform tables:
    some rows carry a JGI numeric id or a free-text gene-name list in the ORF
    column, and those must not be mistaken for locus tags.
    """
    rows, in_tab = [], False
    with open(gpl_path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_tab = True; continue
            if line.startswith("!platform_table_end"):
                break
            if in_tab:
                rows.append(line.rstrip("\n").split("\t"))
    head = rows[0]
    if header_pick is not None:
        orf_col = head.index(header_pick)
    rx = re.compile(pattern) if pattern else None
    out = {}
    for r in rows[1:]:
        if len(r) > orf_col and r[orf_col]:
            v = r[orf_col].strip().strip('"')
            if rx and not rx.match(v):
                continue
            out[r[probe_col]] = v
    return out


def _collapse(expr, probe2gene):
    """Map probes to genes and average probes of the same gene (in log space)."""
    g = pd.Series({p: probe2gene.get(p) for p in expr.index})
    expr = expr.loc[g.notna().values]
    expr = expr.groupby(g.dropna().values).mean()
    expr.index.name = "gene"
    return expr


def _log2(df):
    """Arrays stored as linear intensities -> log2; already-log tables pass through."""
    v = df.to_numpy(dtype=float)
    finite = v[np.isfinite(v)]
    if finite.size and (finite.min() < -0.5 or np.nanmax(finite) < 20):
        return df.astype(float)                     # already log-scale
    return np.log2(df.clip(lower=1).astype(float))


def _cpm_log2(counts):
    lib = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.divide(lib, axis=1) * 1e6 + 1)


# ---------------------------------------------------------------- per dataset
def load_SE7942_ito(geo, arm="WT"):
    tab, meta = _series_matrix(geo / "GSE14225.txt.gz")
    titles = _titles(meta)
    keep, times, series = [], [], []
    pat = re.compile(r"^(WT|kaiABC-null|kaiC-overexpressor) at LL\s*(\d+), rep(\d)$")
    for col, t in zip(tab.columns, titles):
        m = pat.match(t)
        if m and m.group(1) == arm:
            keep.append(col); times.append(float(m.group(2))); series.append("rep" + m.group(3))
    expr = _log2(tab[keep])
    p2g = _probe_map(geo / "GPL8012.annot.txt", header_pick="7942_ID",
                     pattern=r"^Synpcc7942_\w+$")
    expr = _collapse(expr, p2g)
    return expr, pd.DataFrame({"time": times, "series": series}, index=keep)


def load_SE7942_vijayan(geo):
    tab, meta = _series_matrix(geo / "GSE18902.txt.gz")
    titles = _titles(meta)
    keep, times, series = [], [], []
    pat = re.compile(r"^continuous light T = (\d+) hours experiment (\d)$")
    for col, t in zip(tab.columns, titles):
        m = pat.match(t)
        if m:
            keep.append(col); times.append(float(m.group(1))); series.append("exp" + m.group(2))
    expr = _log2(tab[keep])
    p2g = _probe_map(geo / "GPL9534.annot.txt", header_pick="ORF",
                     pattern=r"^Synpcc7942_\w+$")
    expr = _collapse(expr, p2g)
    return expr, pd.DataFrame({"time": times, "series": series}, index=keep)


def load_SYN6803(geo):
    tab, meta = _series_matrix(geo / "GSE79714.txt.gz")
    titles = _titles(meta)
    dayno = {"first": 0, "second": 1, "third": 2}
    pat = re.compile(r"^(\d+) hour after (dark|light) time in the (\w+) day$")
    keep, times = [], []
    for col, t in zip(tab.columns, titles):
        m = pat.match(t)
        if m and m.group(3) in dayno:
            h = float(m.group(1)) + (12.0 if m.group(2) == "dark" else 0.0)
            keep.append(col); times.append(dayno[m.group(3)] * 24 + h)
    expr = _log2(tab[keep])
    p2g = _probe_map(geo / "GPL21668.annot.txt", header_pick="ORF",
                     pattern=r"^(sll|slr|ssl|ssr|sml|smr|ncl|ncr|sll)\d+$")
    expr = _collapse(expr, p2g)
    return expr, pd.DataFrame({"time": times, "series": "s1"}, index=keep)


def load_TERY(geo):
    tab, meta = _series_matrix(geo / "GSE99896.txt.gz")
    titles = _titles(meta)
    zt_line = [l for l in meta["Sample_characteristics_ch1"]
               if l and l[0].startswith("hours into light")][0]
    keep, times, series = [], [], []
    for col, t, z in zip(tab.columns, titles, zt_line):
        if "mixture" in t:
            continue
        m = re.search(r":\s*(\d?)([LD])(\d+)", z)
        if not m:
            continue
        zt = float(m.group(3)) + (12.0 if m.group(2) == "D" else 0.0)
        day = int(m.group(1)) - 1 if m.group(1) else 0   # "2D12" = day 2
        rep = re.search(r"rep(\d)", t)
        keep.append(col); times.append(day * 24 + zt)
        series.append("rep" + (rep.group(1) if rep else "1"))
    expr = _log2(tab[keep])
    p2g = _probe_map(geo / "GPL23571.annot.txt", header_pick="ORF",
                     pattern=r"^Tery_\d+$")
    expr = _collapse(expr, p2g)
    return expr, pd.DataFrame({"time": times, "series": series}, index=keep)


def load_MED4(geo):
    counts = pd.read_csv(geo / "MED4_counts.csv.gz", index_col=0)
    pat = re.compile(r"^(\d+)Hr-(\d+)$")        # drops the -dk-/-lt- perturbation samples
    keep, times = [], []
    for c in counts.columns:
        m = pat.match(c)
        if m:
            keep.append(c); times.append(float(m.group(1)))
    expr = _cpm_log2(counts[keep])
    return expr, pd.DataFrame({"time": times, "series": "s1"}, index=keep)


def load_NPUN(geo):
    counts = pd.read_csv(geo / "NOSTOC_counts.csv.gz", index_col=0)
    counts = counts.drop(columns=[c for c in ("Product", "KEGG1", "KEGG2", "KEGG3")
                                  if c in counts.columns])
    # Dawn/MidLight/Dusk n = day n; MidDark n = the night that follows day n-1,
    # which is where it sits in the sampling order of this series.
    zt = {"Dawn": 0.0, "MidLight": 6.0, "Dusk": 12.0, "MidDark": -6.0}
    pat = re.compile(r"^(Dawn|Dusk|MidLight|MidDark)(\d)-(\d)$")   # bare cols / Dark24 / Light24 dropped
    keep, times, series = [], [], []
    for c in counts.columns:
        m = pat.match(str(c))
        if m:
            keep.append(c); times.append(int(m.group(2)) * 24 + zt[m.group(1)])  # MidDark -> -6 h
            series.append("rep" + m.group(3))
    expr = _cpm_log2(counts[keep].apply(pd.to_numeric, errors="coerce").fillna(0))
    return expr, pd.DataFrame({"time": times, "series": series}, index=keep)


def load_HALO(geo, salt="0.5M"):
    """Halothece PCC 7418 free-running LL. Both biological replicate files are
    used, and the two NaCl concentrations are kept as separate datasets: the
    salinity is an experimental variable of that study, not a batch effect."""
    parts = []
    for rep in ("r1", "r2"):
        t = pd.read_csv(geo / f"HALOTHECE_{rep}.txt.gz", sep="\t", index_col=0)
        t = t[[c for c in t.columns if c.endswith("_RPKM")]]
        t.columns = [f"{c}_{rep}" for c in t.columns]
        parts.append(t)
    tab = pd.concat(parts, axis=1)
    pat = re.compile(r"^([\d.]+M)_LL(\d+)_RPKM_(r\d)$")
    keep, times, series = [], [], []
    for c in tab.columns:
        m = pat.match(c)
        if m and m.group(1) == salt:
            keep.append(c); times.append(float(m.group(2))); series.append(m.group(3))
    expr = np.log2(tab[keep].astype(float) + 1)
    expr.index = [str(i).replace("gene-", "") for i in expr.index]
    return expr, pd.DataFrame({"time": times, "series": series}, index=keep)


LOADERS = {
    "SE7942_ito_WT":       lambda g: load_SE7942_ito(g, "WT"),
    "SE7942_ito_KAINULL":  lambda g: load_SE7942_ito(g, "kaiABC-null"),
    "SE7942_vijayan":      load_SE7942_vijayan,
    "SYN6803":             load_SYN6803,
    "TERY":                load_TERY,
    "MED4":                load_MED4,
    "NPUN":                load_NPUN,
    "HALO_05":             lambda g: load_HALO(g, "0.5M"),
    "HALO_25":             lambda g: load_HALO(g, "2.5M"),
}
