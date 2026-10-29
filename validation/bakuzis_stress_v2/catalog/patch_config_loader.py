"""Patch config_loader.generate_keywords: a per-species BAMAX vector from the variant's own blkdat
(categories/site_index or other) is not a calibrated stand-level cap; skip it with a warning instead of
raising, so the calibrated multipliers still apply. Found 2026-10-29: ACD and IE calibrated arms were
identical to default because the raise aborted the whole keyword block."""
import sys, re
p = sys.argv[1]; t = open(p).read()
old = '''        bamax_values = self._find_bamax_param(cats)
        if bamax_values is not None:
            lines.append(self._format_bamax_keywords(bamax_values, include_comments))
'''
new = '''        bamax_values = self._find_bamax_param(cats)
        if bamax_values is not None:
            try:
                lines.append(self._format_bamax_keywords(bamax_values, include_comments))
            except ValueError as exc:
                # A per-species BAMAX vector (the variant's own blkdat array, e.g. ACD and IE) is not a
                # stand-level cap and is not a calibrated quantity. Skip it and keep the rest of the
                # block; before 2026-10-29 this raise aborted every calibrated keyword for the variant.
                warnings.warn(f"BAMAX skipped for variant {self.variant}: {exc}")
'''
assert t.count(old) == 1, "anchor not found"
t = t.replace(old, new)
if "import warnings" not in t:
    t = t.replace("import numpy as np", "import numpy as np\nimport warnings", 1)
open(p, "w").write(t); print("patched", p)
