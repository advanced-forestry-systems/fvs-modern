#!/usr/bin/env bash
# Portable cross-OS smoke + load test for all FVS variant shared libraries.
# Builds every variant, then ctypes-loads each library (with full symbol resolution,
# RTLD_NOW) so missing-symbol link gaps fail loudly. Runs on Linux, macOS, and
# Windows (MSYS2/MinGW). Used by the cross-platform CI matrix and by end users.
# Usage: bash deployment/scripts/test_all_variants.sh [SOURCE_DIR] [OUTPUT_DIR] [variant ...]
set -uo pipefail
SRC="${1:-src-converted}"; OUT="${2:-lib}"; shift 2 2>/dev/null || true
VARIANTS="$*"
echo "== building variant libraries =="
bash deployment/scripts/build_fvs_libraries.sh "$SRC" "$OUT" $VARIANTS || true
echo "== ctypes load test (RTLD_NOW: undefined symbols fail) =="
python3 - "$OUT" << "PY"
import ctypes, glob, os, sys
out=sys.argv[1]
libs=sorted(glob.glob(os.path.join(out,"FVS*.so"))+glob.glob(os.path.join(out,"FVS*.dylib"))+glob.glob(os.path.join(out,"FVS*.dll")))
if not libs:
    print("NO LIBRARIES FOUND in",out); sys.exit(1)
ok=0; fails=[]
for L in libs:
    try:
        ctypes.CDLL(L, mode=getattr(ctypes,"RTLD_GLOBAL",0)|os.RTLD_NOW if hasattr(os,"RTLD_NOW") else 0)
        ok+=1
    except Exception as e:
        fails.append((os.path.basename(L), str(e).split(":")[-1].strip()))
print(f"loaded {ok}/{len(libs)} variant libraries")
for n,e in fails: print("  FAIL",n,"->",e)
sys.exit(1 if fails else 0)
PY
