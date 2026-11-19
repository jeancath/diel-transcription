"""Gene -> KEGG Orthology mapping, which is the cross-lineage common currency.

Two of the datasets are quantified against current RefSeq locus tags
(TX50_RS..., PCC7418_RS...) while KEGG still keys those genomes on the original
locus tags (PMM..., PCC7418_0001). The RefSeq GFF carries both, so the bridge is
exact rather than name-matched.
"""
import gzip
import re
import urllib.request
from pathlib import Path

KEGG = "https://rest.kegg.jp"


def _cached(url, path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        data = urllib.request.urlopen(url, timeout=120).read()
        path.write_bytes(data)
    return path.read_text(errors="replace")


def gene_to_ko(org, cache_dir):
    """KEGG gene id (without org prefix) -> KO id."""
    txt = _cached(f"{KEGG}/link/ko/{org}", Path(cache_dir) / f"ko_{org}.tsv")
    out = {}
    for line in txt.splitlines():
        if "\t" not in line:
            continue
        g, k = line.split("\t")[:2]
        out[g.split(":", 1)[1]] = k.split(":", 1)[1]
    return out


def ko_names(cache_dir):
    txt = _cached(f"{KEGG}/list/ko", Path(cache_dir) / "ko_names.tsv")
    out = {}
    for line in txt.splitlines():
        if "\t" in line:
            k, v = line.split("\t", 1)
            out[k.replace("ko:", "")] = v
    return out


def locus_tag_bridge(gff_gz):
    """current locus_tag -> old_locus_tag, straight out of the RefSeq GFF."""
    out = {}
    with gzip.open(gff_gz, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            lt = re.search(r"locus_tag=([^;\n]+)", f[8])
            ol = re.search(r"old_locus_tag=([^;\n]+)", f[8])
            if lt and ol:
                out[lt.group(1)] = ol.group(1).split(",")[0]
    return out
