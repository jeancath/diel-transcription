"""Reproducible literature search and screen.

Runs a fixed query set against PubMed, records the number of records retrieved
per query, screens titles and abstracts against stated inclusion criteria, and
writes the whole thing to results/ so the priority statement in the manuscript
can be audited rather than taken on trust.

Inclusion criterion for the priority question: a study that compares diel or
circadian transcriptional programs BETWEEN two or more cyanobacterial lineages
at the level of orthologous genes.
Secondary criterion: a study that addresses the null model used for
cross-species timing comparison.
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from analyze import OUT

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

QUERIES = {
 "Q1 cyano comparative diel": '(cyanobacteri*[tiab]) AND (diel[tiab] OR diurnal[tiab] OR circadian[tiab]) AND (transcriptom*[tiab]) AND (compar*[tiab] OR "across species"[tiab] OR interspecies[tiab] OR conserv*[tiab])',
 "Q2 cyano rhythm meta-analysis": '(cyanobacteri*[tiab]) AND (circadian[tiab] OR diel[tiab]) AND (meta-analysis[tiab] OR reanalysis[tiab] OR "re-analysis"[tiab])',
 "Q3 clock output conservation": '(clock-controlled[tiab] OR "circadian output"[tiab] OR "clock output"[tiab]) AND cyanobacteri*[tiab] AND (conserv*[tiab] OR evolution*[tiab])',
 "Q4 ortholog phase conservation": '(circadian[tiab] OR diel[tiab] OR diurnal[tiab]) AND phase[tiab] AND ortholog*[tiab] AND conserv*[tiab]',
 "Q5 null model timing comparison": '(circadian[tiab] OR diel[tiab] OR diurnal[tiab]) AND ("null model"[tiab] OR permutation[tiab]) AND (phase[tiab] OR timing[tiab]) AND (species[tiab] OR ortholog*[tiab])',
 "Q6 cross-species circadian transcriptome": '(circadian[tiab] OR diurnal[tiab]) AND transcriptom*[tiab] AND ("cross-species"[tiab] OR "between species"[tiab] OR "across species"[tiab])',
 "Q7 rhythm detection reproducibility": '(circadian[tiab] OR rhythmic[tiab]) AND (reproducib*[tiab] OR repeatab*[tiab] OR replicat*[tiab]) AND (transcriptom*[tiab] OR "gene expression"[tiab])',
 "Q8 Synechococcus elongatus circadian transcriptome": '"Synechococcus elongatus"[tiab] AND circadian[tiab] AND (transcriptom*[tiab] OR microarray[tiab] OR "RNA-seq"[tiab])',
}

CYANO = re.compile(r"cyanobact|Synechococcus|Synechocystis|Prochlorococcus|Nostoc|"
                   r"Trichodesmium|Anabaena|Cyanothece|Halothece|Microcystis", re.I)
TIMING = re.compile(r"diel|diurnal|circadian|rhythm|phase|time.?course", re.I)
CROSS = re.compile(r"compar|across species|between species|cross-species|ortholog|"
                   r"lineage|conserv|evolution", re.I)
NULLM = re.compile(r"null model|permutation|background distribution|marginal", re.I)


def get(u, tries=4):
    for i in range(tries):
        try:
            return urllib.request.urlopen(u, timeout=120).read()
        except Exception as e:
            err = e
            time.sleep(3 + 2 * i)
    raise err


def esearch(term, retmax=200):
    d = json.loads(get(E + "esearch.fcgi?" + urllib.parse.urlencode(
        {"db": "pubmed", "term": term, "retmax": retmax, "retmode": "json"})))
    return d["esearchresult"]["count"], d["esearchresult"]["idlist"]


def fetch_records(ids):
    out = {}
    for i in range(0, len(ids), 150):
        txt = get(E + "efetch.fcgi?" + urllib.parse.urlencode(
            {"db": "pubmed", "id": ",".join(ids[i:i + 150]),
             "rettype": "abstract", "retmode": "text"})).decode("utf8", "replace")
        for rec in re.split(r"\n\n(?=\d+\. )", txt):
            m = re.search(r"PMID:\s*(\d+)", rec)
            if m:
                out[m.group(1)] = rec
        time.sleep(0.4)
    return out


FERRARI_PMID = "30760717"   # Ferrari et al. 2019, Nat Commun 10:737


def citing_screen():
    """Screen every record that cites the closest prior study.

    One-level citation traversal via PubMed `elink` (pubmed_pubmed_citedin).
    Regenerated from the live citation graph, so the count can grow after the
    search date recorded in LITERATURE-SEARCH.md; that date, not this run, is
    what the manuscript cites.
    """
    d = json.loads(get(E + "elink.fcgi?" + urllib.parse.urlencode(
        {"dbfrom": "pubmed", "db": "pubmed", "id": FERRARI_PMID,
         "linkname": "pubmed_pubmed_citedin", "retmode": "json"})))
    ids = []
    for ls in d["linksets"]:
        for db in ls.get("linksetdbs", []):
            if db.get("linkname") == "pubmed_pubmed_citedin":
                ids = [str(x) for x in db["links"]]
    print(f"citing records for PMID {FERRARI_PMID}: {len(ids)}", flush=True)
    recs = fetch_records(sorted(ids))
    rows = []
    for pmid, rec in recs.items():
        ti = re.search(r"\n\n([^\n].*?)\n\n", rec, re.S)
        title = " ".join(ti.group(1).split()) if ti else ""
        rows.append(dict(pmid=pmid, title=title[:140],
                         cyano=bool(CYANO.search(rec)),
                         comparative=bool(CROSS.search(rec)),
                         timing=bool(TIMING.search(rec)),
                         nullmodel=bool(NULLM.search(rec))))
    t = pd.DataFrame(rows)
    t["flag_fulltext"] = t.cyano & t.timing
    t.to_csv(OUT / "lit_ferrari_citing_screen.csv", index=False)
    print(f"  cyanobacterial: {int(t.cyano.sum())}   null-model: "
          f"{int(t.nullmodel.sum())}   flagged: {int(t.flag_fulltext.sum())}")
    return t


def main():
    today = date.today().isoformat()
    counts, all_ids = [], set()
    for name, term in QUERIES.items():
        c, ids = esearch(term)
        counts.append(dict(query_id=name, query=term, date=today,
                           records_retrieved=int(c), records_screened=len(ids)))
        all_ids |= set(ids)
        print(f"{name:42s} {c:>5} records", flush=True)
        time.sleep(0.4)
    pd.DataFrame(counts).to_csv(OUT / "lit_query_log.csv", index=False)

    recs = fetch_records(sorted(all_ids))
    rows = []
    for pmid, rec in recs.items():
        ti = re.search(r"\n\n([^\n].*?)\n\n", rec, re.S)
        title = " ".join(ti.group(1).split()) if ti else ""
        rows.append(dict(pmid=pmid, title=title[:140],
                         cyano=bool(CYANO.search(rec)), timing=bool(TIMING.search(rec)),
                         cross=bool(CROSS.search(rec)), nullmodel=bool(NULLM.search(rec))))
    t = pd.DataFrame(rows)
    t["meets_primary"] = t.cyano & t.timing & t.cross
    t["meets_secondary"] = t.timing & t.cross & t.nullmodel
    t.to_csv(OUT / "lit_screen.csv", index=False)
    print(f"\nunique records screened: {len(t)}")
    print(f"meet primary criterion (cyanobacteria + timing + cross-lineage): "
          f"{int(t.meets_primary.sum())}")
    print(f"meet secondary criterion (timing + cross-species + null model): "
          f"{int(t.meets_secondary.sum())}")
    if t.meets_primary.any():
        print("\n-- flagged for full text (primary):")
        print(t[t.meets_primary][["pmid", "title"]].to_string(index=False))
    if t.meets_secondary.any():
        print("\n-- flagged for full text (secondary):")
        print(t[t.meets_secondary][["pmid", "title"]].to_string(index=False))

    print("\n== citation network of the closest prior study")
    citing_screen()


if __name__ == "__main__":
    main()
