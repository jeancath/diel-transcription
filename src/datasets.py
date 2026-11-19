"""Dataset registry for the cross-lineage diel/circadian transcriptome meta-analysis.

Each entry says where the expression data lives, how to recover the sampling time
(hours) for every column, and which KEGG organism the gene IDs belong to.
Nothing here is inferred at runtime: sampling times come from the GEO sample
titles / characteristics, which is the whole reason this panel is usable.
"""
from dataclasses import dataclass, field

GEO_DIR = None  # set by driver


@dataclass
class Dataset:
    did: str
    gse: str
    organism: str          # display name
    lineage: str           # short label used in figures
    kegg: str              # KEGG organism code
    light: str             # "LL" = free-running constant light, "LD" = driven light/dark
    platform: str          # "array" or "rnaseq"
    source: str            # file in the geo dir
    gpl: str = ""          # platform annotation file (arrays)
    note: str = ""
    role: str = "test"     # "test", "control" (clock-null), "replicate_study"


DATASETS = [
    # --- Synechococcus elongatus PCC 7942: the model strain, three independent studies
    Dataset("SE7942_ito_WT", "GSE14225", "Synechococcus elongatus PCC 7942",
            "S. elongatus", "syf", "LL", "array", "GSE14225.txt.gz", "GPL8012.annot.txt",
            "Ito et al. 2009; WT, LL4-LL52 every 4 h, 2 biological series"),
    Dataset("SE7942_ito_KAINULL", "GSE14225", "Synechococcus elongatus PCC 7942 kaiABC-null",
            "S. elongatus (kai-null)", "syf", "LL", "array", "GSE14225.txt.gz", "GPL8012.annot.txt",
            "Ito et al. 2009; kaiABC-null arm - negative control for the rhythm caller",
            role="control"),
    Dataset("SE7942_vijayan", "GSE18902", "Synechococcus elongatus PCC 7942",
            "S. elongatus (study 2)", "syf", "LL", "array", "GSE18902.txt.gz", "GPL9534.annot.txt",
            "Vijayan et al. 2009; LL every 4 h, 2 experiments, different platform",
            role="replicate_study"),
    # --- other lineages
    Dataset("SYN6803", "GSE79714", "Synechocystis sp. PCC 6803",
            "Synechocystis", "syn", "LD", "array", "GSE79714.txt.gz", "GPL21668.annot.txt",
            "LD 12:12, every 2 h over 2 days"),
    Dataset("TERY", "GSE99896", "Trichodesmium erythraeum IMS101",
            "Trichodesmium", "ter", "LD", "array", "GSE99896.txt.gz", "GPL23571.annot.txt",
            "LD, every 3 h, 3 replicates"),
    Dataset("MED4", "GSE314951", "Prochlorococcus marinus MED4",
            "Prochlorococcus", "pmm", "LD", "rnaseq", "MED4_counts.csv.gz",
            note="LD, every 2 h; kaiA-less hourglass lineage"),
    Dataset("NPUN", "GSE275682", "Nostoc punctiforme PCC 73102",
            "Nostoc", "npu", "LD", "rnaseq", "NOSTOC_counts.csv.gz",
            note="diel, 4 phases/day over 3 days, 3 replicates"),
    Dataset("HALO_05", "GSE205640", "Halothece sp. PCC 7418 (0.5 M NaCl)",
            "Halothece 0.5M", "hao", "LL", "rnaseq", "HALOTHECE_r*.txt.gz",
            note="free-running LL every 6 h over 30 h, 2 replicates, low salt"),
    Dataset("HALO_25", "GSE205640", "Halothece sp. PCC 7418 (2.5 M NaCl)",
            "Halothece 2.5M", "hao", "LL", "rnaseq", "HALOTHECE_r*.txt.gz",
            note="free-running LL every 6 h over 30 h, 2 replicates, high salt"),
]

BY_ID = {d.did: d for d in DATASETS}
