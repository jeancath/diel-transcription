"""Where downloaded inputs live.

Set DIEL_DATA_DIR to put the ~25 MB of downloaded GEO/NCBI/KEGG files wherever
you like; it defaults to a `data/` directory beside this repository, so the
pipeline runs on any machine without editing code.
"""
import os
from pathlib import Path

DATA = Path(os.environ.get(
    "DIEL_DATA_DIR",
    Path(__file__).resolve().parents[1] / "data")).expanduser()
GEO = DATA / "geo"
KEGG = DATA / "kegg"
RBH = DATA / "rbh"
PROTEOMES = DATA / "proteomes"

for _d in (GEO, KEGG, RBH, PROTEOMES):
    _d.mkdir(parents=True, exist_ok=True)
