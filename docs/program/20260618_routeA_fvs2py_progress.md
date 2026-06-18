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

## Rebuild result (2026-06-18, this session)

Rebuilt FVSne fresh into lib-test/FVSne.so from the current source with the build-fixes build script
(env FC=gfortran CC=gcc build_fvs_libraries.sh src-converted ./lib-test ne). Results, verified in process:

- The keyrdr.f90 EOF blocker is GONE. The rebuilt .so reads an in-process keyword file without the fort.15
  EOF error.
- In-process load, load_keyfile, and run to stop point 7 work: the stand loads and dims['ntrees'] returns
  the correct live tree count (67 for the test stand).
- The new get_tree_attr and set_tree_attr bindings WORK in process: reading 'dbh' and 'dg' returns
  per-tree arrays of the right length, and set_tree_attr('dg', ...) executes without error.

The one remaining blocker, now isolated to the engine itself (not the wrapper, not the tree-attr code, not
the database): stepping the in-process engine to stop point 5 (after growth and mortality are computed,
before applied) SEGFAULTS inside the FVS step routine (_base.py run, the self._fvs call). This was
reproduced four ways, which together exclude the obvious causes:

- Trimmed input-only keyword file: segfault at the first stop-point-5 step.
- Full keyword file with DSNOUT on a separate output database (so the output DB opens cleanly): still
  segfaults at the first stop-point-5 step. This rules out the database.
- Without touching any tree attribute (no get/set before the step): still segfaults. This rules out the
  new tree-attr binding.
- At stop point 7 the tree count is correct (67) but per-tree dbh and dg read as 0.0, confirming the
  working tree arrays are not yet populated that early; the crash is at the growth step regardless.

The earlier "unrecognized token: '" output-DB error appeared only when DSNIN and DSNOUT shared one file;
pointing DSNOUT at a separate file clears it, after which the run reaches growth and then segfaults.

Conclusion: the in-process growth / stop-point-5 path in the .so faults. This is a source-level engine bug
(the stop-point restart wiring or the growth initialization under the in-process API), not a Python-side
or data issue. Fixing it needs gdb on the loaded .so to capture the faulting Fortran frame.

Next steps to finish Route A (revised):

1. Run FVSne under gdb in process: load the .so, set a breakpoint, step to stop point 5, and capture the
   faulting frame and backtrace. Likely candidates are the stop-point restart save/restore (the SPESET /
   restart-code machinery) or a growth/FFE initializer that the standalone executable reaches differently.
2. Compare the in-process stop-point path against rFVS, which exercises the same fvs API stop points in R;
   if rFVS steps stop point 5 cleanly on the same .so, the gap is in how fvs2py drives the run loop, not
   the engine.
3. If the stop-point path proves unreliable, fall back to a cycle-boundary override: run a full cycle, read
   the projected tree list, rescale diameters to the fvs-conus prediction, reload, and project the next
   cycle. Coarser than a mid-cycle override but avoids the stop-point-5 fault entirely.

Earlier framing of the remaining work (kept for reference):

## Next steps to finish Route A

1. Read and modify tree attrs only after the stand is fully set up. Step past stop point 7 to the first
   point where the tree arrays are populated (run one full cycle with no stop, or stop at code 1/2 at the
   first cycle), confirm get_tree_attr('dbh') returns the real diameters, then introduce the stop-point-5
   override. The zero-valued reads at stop point 7 are the cause of the growth-step segfault.
2. Resolve the output-DB path for the full template: point DSNOUT at a separate file from DSNIN, or drop
   the output-DB keywords and read all results through the summary and tree_table APIs (the latter is the
   intended in-process pattern). This clears the "unrecognized token: '" output-DB error.
3. If DB-driven loading stays fragile, drive the engine entirely through the API: fvsAddTrees to load the
   stand in memory (exported and bound), bypassing the keyword-file DSN altogether. Cleaner long-term path.
4. Once a stand projects cleanly in process, wire the fvs-conus DG and HD predictions into the
   stop-point-5 loop (replace the +30 percent stand-in with the fvs-conus equation evaluation), producing
   true in-engine arms C and D, and drop the emulation caveat from the four-arm.

## Artifacts

- deployment/fvs2py/fvs2py/_base.py: get_tree_attr / set_tree_attr / tree_table (the unblock).
- diagnostics_2026-06-16/test_treeattr_injection.py: the in-process injection test (loads a stand, stops
  at code 5, overwrites dg, resumes). Currently blocked at the keyrdr EOF described above.

Bottom line: the maintainer-level capability gap (in-process per-tree read and write) is closed and
committed. What remains is a build or stand-loading-path fix so a stand projects in process, after which
the true in-engine four-arm is a short step. Until then, the committed four-arm uses the per-species
multiplier emulation, which is honestly labeled.
