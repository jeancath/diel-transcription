"""Record exactly what produced the current results, and verify it later.

Writes results/PROVENANCE.json: the commit, whether the tree was clean, the
interpreter and package versions, the platform, a SHA-256 for every raw input
and every generated output, and the declared producing target for each result
file.

`--verify` re-reads the manifest and reports any input whose checksum has
changed since, any output missing, and any output not covered by the manifest.
An input that changed invalidates the outputs derived from it, which a build
log alone would not reveal.

Run: python3 src/provenance.py            (write)
     python3 src/provenance.py --verify   (check)
"""
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from paths import DATA
from manuscript_numbers import _PRODUCER

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "figures"
MANIFEST = RES / "PROVENANCE.json"
INPUT_DIRS = ("geo", "kegg", "proteomes", "rbh")


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def git(*args, default="unknown"):
    try:
        r = subprocess.run(["git", *args], cwd=ROOT.parent,
                           capture_output=True, text=True, timeout=30)
        return r.stdout.strip() if r.returncode == 0 else default
    except Exception:
        return default


def producer_of(name):
    for k, t in _PRODUCER.items():
        if name.startswith(k):
            return f"make {t}"
    return None


def inputs():
    out = {}
    for sub in INPUT_DIRS:
        d = DATA / sub
        if not d.exists():
            continue
        for f in sorted(d.rglob("*")):
            if f.is_file():
                out[f"{sub}/{f.relative_to(d)}"] = {
                    "sha256": sha256(f), "bytes": f.stat().st_size}
    return out


def outputs():
    out = {}
    for f in sorted(list(RES.glob("*.csv")) + list(RES.glob("*.csv.gz"))
                    + list(RES.glob("*.json")) + list(FIG.glob("*.png"))):
        if f.name in ("PROVENANCE.json",) or f.name.endswith(".partial.csv"):
            continue
        rel = f"{f.parent.name}/{f.name}"
        out[rel] = {"sha256": sha256(f), "bytes": f.stat().st_size,
                    "producer": producer_of(f.name)}
    return out


def versions():
    v = {"python": sys.version.split()[0], "platform": platform.platform(),
         "machine": platform.machine()}
    for mod in ("numpy", "pandas", "scipy", "matplotlib", "Bio", "pyswrd"):
        try:
            m = __import__(mod)
            v[mod] = getattr(m, "__version__", "unknown")
        except Exception as e:
            v[mod] = f"unavailable: {e}"
    return v


def write():
    man = {
        "written_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git": {
            "commit": git("rev-parse", "HEAD"),
            "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
            "tree_clean": git("status", "--porcelain") == "",
            "dirty_paths": [l[3:] for l in git("status", "--porcelain").splitlines()][:40],
            "describe": git("describe", "--always", "--dirty"),
        },
        "versions": versions(),
        "data_dir": str(DATA),
        "inputs": inputs(),
        "outputs": outputs(),
    }
    man["summary"] = {
        "n_inputs": len(man["inputs"]),
        "input_bytes": sum(d["bytes"] for d in man["inputs"].values()),
        "n_outputs": len(man["outputs"]),
        "outputs_without_producer": sorted(
            k for k, d in man["outputs"].items() if d["producer"] is None),
    }
    MANIFEST.write_text(json.dumps(man, indent=2, sort_keys=True) + "\n")
    s = man["summary"]
    print(f"wrote {MANIFEST.relative_to(ROOT)}")
    print(f"  commit      {man['git']['describe']}  (clean tree: {man['git']['tree_clean']})")
    print(f"  inputs      {s['n_inputs']} files, {s['input_bytes'] / 1e6:.1f} MB")
    print(f"  outputs     {s['n_outputs']} files")
    if s["outputs_without_producer"]:
        print(f"  NO PRODUCER {s['outputs_without_producer']}")
    return man


def verify():
    if not MANIFEST.exists():
        print("no PROVENANCE.json; run `python3 src/provenance.py` first")
        return 1
    man = json.loads(MANIFEST.read_text())
    problems = []
    now_in = inputs()
    for k, d in man["inputs"].items():
        if k not in now_in:
            problems.append(f"input missing: {k}")
        elif now_in[k]["sha256"] != d["sha256"]:
            problems.append(f"input CHANGED since the manifest: {k}")
    for k in sorted(set(now_in) - set(man["inputs"])):
        problems.append(f"input present but not in the manifest: {k}")
    now_out = outputs()
    for k, d in man["outputs"].items():
        if k not in now_out:
            problems.append(f"output missing: {k}")
        elif now_out[k]["sha256"] != d["sha256"]:
            problems.append(f"output changed since the manifest: {k}")
    for k in sorted(set(now_out) - set(man["outputs"])):
        problems.append(f"output present but not in the manifest: {k}")
    missing_prod = [k for k, d in man["outputs"].items() if d["producer"] is None]
    for k in missing_prod:
        problems.append(f"output has no declared producer: {k}")
    print(f"manifest written {man['written_utc']} at commit "
          f"{man['git']['describe']}")
    print(f"  {len(man['inputs'])} inputs, {len(man['outputs'])} outputs checked")
    if problems:
        print(f"\nPROVENANCE MISMATCH: {len(problems)} problem(s)")
        for p in problems[:40]:
            print(f"  - {p}")
        return 1
    print("\nPROVENANCE OK: every input and output matches the manifest")
    return 0


if __name__ == "__main__":
    sys.exit(verify() if "--verify" in sys.argv else (write() and 0))
