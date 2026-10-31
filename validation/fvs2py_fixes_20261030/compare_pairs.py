"""Compare CSV pairs: usage compare_pairs.py out.csv label a.csv b.csv [label a b ...]"""
import sys
import numpy as np, pandas as pd
out = []
args = sys.argv[2:]
for i in range(0, len(args), 3):
    lab, fa, fb = args[i:i + 3]
    a, b = pd.read_csv(fa), pd.read_csv(fb)
    common = [c for c in a.columns if c in b.columns]
    num = [c for c in common if np.issubdtype(a[c].dtype, np.number) and np.issubdtype(b[c].dtype, np.number)]
    same_shape = len(a) == len(b)
    mad = float(np.nanmax(np.abs(a[num].to_numpy(float) - b[num].to_numpy(float)))) if same_shape and num else float("nan")
    out.append(dict(comparison=lab, file_a=fa.split("/")[-1], file_b=fb.split("/")[-1], rows_a=len(a), rows_b=len(b),
                    cols_compared=len(common), identical=bool(same_shape and a[common].equals(b[common])),
                    max_abs_diff=mad, only_in_a=";".join(c for c in a.columns if c not in b.columns),
                    only_in_b=";".join(c for c in b.columns if c not in a.columns)))
d = pd.DataFrame(out); d.to_csv(sys.argv[1], index=False); print(d.to_string(index=False))
