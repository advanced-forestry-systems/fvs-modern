# Route A (fvs2py in-engine injection): progress and the next obstacle
2026-06-18. Approved by the PI. Goal: run the fvs-conus equations inside the FVS engine so arms C and D
are a true single-framework comparison rather than the per-species BAIMULT emulation.

## What was the actual blocker, and what is now unblocked

The documented blocker, "fvs2py in-process tree loading," was concretely this: fvs2py wrapped species
attributes, stop points, and the run loop, but it never wrapped fvsTreeAttr, so it could not read or write
per-tree state (dbh, dg, ht, htg) in the engine's memory. Without that, there was no way to overwrite the
diameter growth the engine computed with an external prediction.

Removed this session. fvs2py now exposes per-tree access (deployment/fvs2py/fvs2py/_base.py):

- get_tree_attr(attr): read a per-tree attribute (length = live tree count, read fresh each call)
- set_tree_attr(attr, arr): write a per-tree attribute
- tree_table(): read several attributes into a DataFrame

These wrap the fvsTreeAttr C entry point (already exported by every variant .so and already bound as
self._fvsTreeAttr in _core.py). The injection mechanism is now complete in principle: run to FVS restart
code 5 (after diameter growth and mortality are computed but before they are applied), read dg, overwrite
with the fvs-conus prediction, write it back with set_tree_attr, and resume so the engine applies it.

Verified: FVSne.so loads in process via ctypes and the stop-point mechanism (code 5, year -1) is reachable.

## The next obstacle (precise)

In-process keyword-file loading hits a Fortran runtime error in base/keyrdr.f90 line 47, "Sequential READ
or WRITE not allowed after EOF marker," on unit 15 (fort.15). This is the same signature as the historical
PN/SN/IE keyrdr EOF blocker recorded in KNOWN_ISSUES, which was fixed on branch build-fixes-2026-05-06 for
the standalone executables. The standalone FVS executable reads the same keyword file fine via subprocess;
the in-process .so path through fvsSetCmdLine(--keywordfile=...) trips the EOF. The most likely cause is
that the installed FVSne.so was not rebuilt with the keyrdr fix, or the in-process keyword-reader needs the
stub handling that the build-fixes branch added.

## Next steps to finish Route A

1. Rebuild the variant .so files from the build-fixes-2026-05-06 source (the keyrdr / INCLUDE-file and stub
   fixes) so the in-process keyword reader does not hit the EOF, and confirm with the tree-attr test.
2. Alternatively, drive the engine entirely through the API (fvsAddTrees to load the stand in memory rather
   than through a keyword-file DSN), which sidesteps the keyword-file reader. fvsAddTrees is exported and
   bound; this is the cleaner long-term path for a pure in-process harness.
3. Once a stand projects in process, wire the fvs-conus DG and HD predictions into the stop-point-5 loop
   (replace the +30 percent stand-in in test_treeattr_injection.py with the fvs-conus equation evaluation),
   producing true in-engine arms C and D, and drop the emulation caveat from the four-arm.

## Artifacts

- deployment/fvs2py/fvs2py/_base.py: get_tree_attr / set_tree_attr / tree_table (the unblock).
- diagnostics_2026-06-16/test_treeattr_injection.py: the in-process injection test (loads a stand, stops
  at code 5, overwrites dg, resumes). Currently blocked at the keyrdr EOF described above.

Bottom line: the maintainer-level capability gap (in-process per-tree read and write) is closed and
committed. What remains is a build or stand-loading-path fix so a stand projects in process, after which
the true in-engine four-arm is a short step. Until then, the committed four-arm uses the per-species
multiplier emulation, which is honestly labeled.
