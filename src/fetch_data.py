"""Download every input this analysis uses. All public, no credentials.

GEO series matrices carry the sample metadata (and, for the microarray studies,
the expression table); the RNA-seq studies keep their count tables as GEO
supplementary files. The two NCBI GFFs supply the current->original locus tag
bridge that KEGG still keys those genomes on.
"""
import sys
import urllib.request
from pathlib import Path

from paths import GEO
FTP = "https://ftp.ncbi.nlm.nih.gov"

SERIES = ["GSE14225", "GSE18902", "GSE79714", "GSE99896",
          "GSE314951", "GSE275682", "GSE205640"]
SUPPL = {
    "MED4_counts.csv.gz":   "/geo/series/GSE314nnn/GSE314951/suppl/GSE314951_rnaseq_counts.csv.gz",
    "NOSTOC_counts.csv.gz": "/geo/series/GSE275nnn/GSE275682/suppl/GSE275682_sansTRNA.csv.gz",
    "HALOTHECE_r1.txt.gz":  "/geo/series/GSE205nnn/GSE205640/suppl/GSE205640_0.5M_2.5M_rep1.txt.gz",
    "HALOTHECE_r2.txt.gz":  "/geo/series/GSE205nnn/GSE205640/suppl/GSE205640_0.5M_2.5M_rep2.txt.gz",
}
GFF = {
    "MED4.gff.gz": "/genomes/all/GCF/000/011/465/GCF_000011465.1_ASM1146v1/"
                   "GCF_000011465.1_ASM1146v1_genomic.gff.gz",
    "HAO.gff.gz":  "/genomes/all/GCF/000/317/635/GCF_000317635.1_ASM31763v1/"
                   "GCF_000317635.1_ASM31763v1_genomic.gff.gz",
}
PLATFORMS = ["GPL8012", "GPL9534", "GPL21668", "GPL23571"]


def grab(url, dest):
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return
    print("fetch", dest.name, flush=True)
    dest.write_bytes(urllib.request.urlopen(url, timeout=300).read())


# ---- genome annotation and proteomes for the RBH orthology arm
#
# These five genomes' protein FASTA, feature table and GFF were originally
# downloaded by hand, so a clean checkout could not reproduce the RBH stage or
# the independence check - `make all` failed on a missing data/rbh/SE7942.gff.gz.
# The accessions are read from the #!genome-build-accession line of the GFFs
# that were in use, so this fetches the same assemblies, not newer ones.
RBH_ASSEMBLIES = {
    "SE7942":  ("GCF_000012525.1", "ASM1252v1"),   # Synechococcus elongatus PCC 7942
    "SYN6803": ("GCF_000009725.1", "ASM972v1"),    # Synechocystis sp. PCC 6803
    "NPUN":    ("GCF_000020025.1", "ASM2002v1"),   # Nostoc punctiforme PCC 73102
    "TERY":    ("GCF_000014265.1", "ASM1426v1"),   # Trichodesmium erythraeum IMS101
    "MED4":    ("GCF_000011465.1", "ASM1146v1"),   # Prochlorococcus marinus MED4
}
RBH_SUFFIXES = {"faa.gz": "protein.faa.gz",
                "ft.gz": "feature_table.txt.gz",
                "gff.gz": "genomic.gff.gz"}


def _assembly_dir(acc):
    """NCBI's six-digit-triplet layout, e.g. GCF_000012525.1 -> 000/012/525."""
    digits = acc.split("_")[1].split(".")[0]
    return "/".join(digits[i:i + 3] for i in range(0, 9, 3))


def fetch_rbh_inputs():
    from paths import RBH
    RBH.mkdir(parents=True, exist_ok=True)
    for key, (acc, asm) in RBH_ASSEMBLIES.items():
        base = (f"/genomes/all/{acc.split('_')[0]}/{_assembly_dir(acc)}/"
                f"{acc}_{asm}/{acc}_{asm}_")
        for local, remote in RBH_SUFFIXES.items():
            dest = RBH / f"{key}.{local}"
            if dest.exists() and dest.stat().st_size > 0:
                continue
            print(f"fetch {dest.name}", flush=True)
            grab(FTP + base + remote, dest)


def main():
    GEO.mkdir(parents=True, exist_ok=True)
    for g in SERIES:
        pre = g[:-3] + "nnn"
        grab(f"{FTP}/geo/series/{pre}/{g}/matrix/{g}_series_matrix.txt.gz", GEO / f"{g}.txt.gz")
    for name, path in {**SUPPL, **GFF}.items():
        grab(FTP + path, GEO / name)
    for p in PLATFORMS:
        grab(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={p}&targ=self&form=text&view=full",
             GEO / f"{p}.annot.txt")
    fetch_rbh_inputs()
    print("all inputs present in", GEO.parent)


if __name__ == "__main__":
    sys.exit(main())
