"""Checkpointing and per-cell seeding for the long simulations.

Two problems this solves.

1. A simulation that writes its CSV only at the end loses everything if the
   process is interrupted. Each cell of the design grid is appended to a
   `.partial.csv` as soon as it finishes, so a restart resumes.

2. A single RNG stream advanced across cells makes a resumed run differ from
   an uninterrupted one, because the skipped cells no longer consume draws.
   `cell_rng` derives an independent stream from the seed and the cell's own
   coordinates, so a cell's result depends only on its coordinates - the same
   whether it is computed first, last, or after a restart.
"""
import json
import os
import platform
import subprocess
import sys
import tempfile
import time
import numpy as np
import pandas as pd
from pathlib import Path


def _git_commit():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"],
                              cwd=Path(__file__).resolve().parents[2],
                              capture_output=True, text=True, timeout=20).stdout.strip()
    except Exception:
        return "unknown"


def cell_rng(seed, key):
    """Independent generator for one grid cell.

    `key` is the cell's coordinates. Floats are scaled to integers so that
    e.g. 0.4 and 0.6 give different streams; strings are hashed stably.
    """
    parts = [int(seed)]
    for k in key:
        if isinstance(k, str):
            parts.append(int.from_bytes(k.encode()[:8].ljust(8, b"\0"), "big") % (2 ** 31))
        else:
            parts.append(int(round(float(k) * 1000)))
    return np.random.default_rng(parts)


class Checkpoint:
    """Append-per-cell CSV with resume."""

    def __init__(self, out_dir, final_name, key_cols):
        self.final = Path(out_dir) / final_name
        self.partial = Path(out_dir) / (final_name.replace(".csv", "") + ".partial.csv")
        self.key_cols = list(key_cols)
        self.rows = []
        self._failures = []
        self._done = set()
        self._t0 = time.time()
        self.failures_path = Path(out_dir) / (final_name.replace(".csv", "") + ".failures.csv")
        self.manifest_path = Path(out_dir) / (final_name.replace(".csv", "") + ".manifest.json")
        if self.partial.exists():
            prev = pd.read_csv(self.partial)
            self.rows = prev.to_dict("records")
            self._done = {self._key(r) for r in self.rows}

    def _key(self, row):
        return tuple(str(row[c]) for c in self.key_cols)

    def done(self, key):
        return tuple(str(k) for k in key) in self._done

    def _write_atomic(self, path, text):
        """Write via a temporary file in the same directory, then rename.

        A direct write that is killed part-way leaves a truncated CSV that
        still parses - a short final row looks like a valid result. os.replace
        is atomic on POSIX, so a reader sees either the previous file or the
        complete new one, never a partial one.
        """
        d = path.parent
        fd, tmp = tempfile.mkstemp(dir=d, prefix=path.name + ".", suffix=".tmp")
        try:
            with os.fdopen(fd, "w") as f:
                f.write(text)
                f.flush()
                os.fsync(f.fileno())
            os.replace(tmp, path)
        except BaseException:
            Path(tmp).unlink(missing_ok=True)
            raise

    def add(self, row):
        self.rows.append(row)
        self._done.add(self._key(row))
        self._write_atomic(self.partial, pd.DataFrame(self.rows).to_csv(index=False))
        print(row, flush=True)

    def failed(self, key, error):
        """Record a cell that raised, so a gap in the grid is explicit."""
        self._failures.append({**dict(zip(self.key_cols, key)), "error": str(error)[:300]})
        self._write_atomic(self.failures_path,
                           pd.DataFrame(self._failures).to_csv(index=False))
        print(f"CELL FAILED {key}: {error}", flush=True)

    def finish(self, sort_by=None, expected=None, config=None):
        """Write the final CSV plus a manifest of what ran.

        `expected` is the full set of cell keys the grid was meant to cover.
        Anything missing from the output is reported as skipped rather than
        left as a silent gap.
        """
        t = pd.DataFrame(self.rows)
        if sort_by:
            t = t.sort_values(sort_by).reset_index(drop=True)
        self._write_atomic(self.final, t.to_csv(index=False))
        self.partial.unlink(missing_ok=True)
        done = {self._key(r) for r in self.rows}
        skipped = []
        if expected is not None:
            skipped = [list(k) for k in expected
                       if tuple(str(x) for x in k) not in done]
        man = {
            "output": self.final.name,
            "script": Path(sys.argv[0]).name,
            "cells_completed": len(self.rows),
            "cells_failed": len(self._failures),
            "cells_skipped": skipped,
            "cells_expected": (len(expected) if expected is not None else None),
            "key_columns": self.key_cols,
            "config": config or {},
            "git_commit": _git_commit(),
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "platform": platform.platform(),
            "wall_seconds": round(time.time() - self._t0, 1),
        }
        self._write_atomic(self.manifest_path, json.dumps(man, indent=2) + "\n")
        print(f"\nwrote {self.final.name} "
              f"({len(self.rows)} cells, {len(self._failures)} failed, "
              f"{len(skipped)} skipped); manifest in {self.manifest_path.name}",
              flush=True)
        return t
