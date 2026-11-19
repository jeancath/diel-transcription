"""Render MANUSCRIPT.md from MANUSCRIPT.template.md + results/.

Every statistic in the manuscript is a {{placeholder}} resolved from
src/manuscript_numbers.py, which reads the result CSVs. Tables marked
{{table:name}} are generated whole. Running this after any analysis change
updates the manuscript; src/check_manuscript.py fails if the two have drifted.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import manuscript_numbers as MN

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "MANUSCRIPT.template.md"
TARGET = ROOT / "MANUSCRIPT.md"
# Case-sensitive keys: kaiB_amp, stratA_sig, drop_MED4. A lowercase-only
# character class here silently leaves such a placeholder as literal text in
# the output - it is not substituted and it is not reported missing, because
# the regex never recognises it. check_manuscript.py therefore also scans the
# rendered file for residual braces, independently of this pattern.
PLACEHOLDER = re.compile(r"\{\{([A-Za-z0-9_:]+)\}\}")


def render(text=None, values=None, tables=None):
    """Render the template. `values`/`tables` are injectable so that
    src/audit_renderer.py can test which keys actually reach the output."""
    text = text if text is not None else TEMPLATE.read_text()
    values = MN.collect() if values is None else values
    tables = MN.tables(values) if tables is None else tables
    missing = []

    def sub(m):
        key = m.group(1)
        if key.startswith("table:"):
            name = key.split(":", 1)[1]
            if name not in tables:
                missing.append(key)
                return m.group(0)
            return tables[name]
        if key not in values:
            missing.append(key)
            return m.group(0)
        return str(values[key])

    out = PLACEHOLDER.sub(sub, text)
    if missing:
        raise KeyError(f"unresolved placeholders: {sorted(set(missing))}")
    # independent of PLACEHOLDER: nothing brace-delimited may survive
    stray = re.findall(r"\{\{[^}]{0,80}\}\}", out)
    if stray:
        raise KeyError(f"residual template braces in output: {sorted(set(stray))}")
    return out


if __name__ == "__main__":
    TARGET.write_text(render())
    print(f"rendered {TARGET.name} from {TEMPLATE.name}")
