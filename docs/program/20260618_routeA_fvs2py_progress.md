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

Conclusion (after gdb): the fault is a data-loading chain, not a growth-equation bug. gdb on the loaded
.so captured the faulting frame:

    #0 extree () at base/extree.f90:38   ->  IMCI = IMC(INS1)   (INS1 = INS(I), I = 1..6)
    #1 fvs (irtncd=0) at base/fvs.f90:306  ->  CALL EXTREE ("assign the example trees to the output arrays")

extree assigns up to six example trees to the output arrays every cycle; it indexes IMC and ISP by the
example-tree indices INS(1..6). In process those indices are unset (0), so IMC(0) reads out of bounds and
segfaults. INS is set during the cycle's tree distribution only after the stand's trees are fully loaded.
The trees never load because the in-process database read fails first: the recurring SQLite error
"unrecognized token: '" means the in-process DBS extension cannot execute the STANDSQL / TREESQL against
the SQLite file (the error persists even with the literal stand id substituted for the %StandID% macro and
with TREELIST removed). So the true root cause is the in-process SQLite read; the extree segfault is a
downstream symptom of an unpopulated stand.

Full chain: in-process DBS SQL fails (unrecognized token) -> trees not loaded -> example-tree indices INS
unset -> extree.f90:38 IMC(INS(I)) out of bounds -> SIGSEGV.

Next steps to finish Route A (revised after the gdb root-cause):

1. RECOMMENDED: load the stand through the API with fvsAddTrees instead of the keyword-file database DSN.
   The blocker is the in-process DBS SQLite read, so bypassing the database removes the whole failure
   chain. fvsAddTrees is exported and bound in _core.py; its signature (base/apisubs.f90:844), to wrap as
   a fvs2py method, is:

       fvsAddTrees(in_dbh, in_species, in_ht, in_cratio, in_plot, in_tpa, ntrees, rtnCode)

   all six inputs real(kind=8) arrays of length ntrees passed by reference; the routine appends them to
   the engine arrays (dbh, isp, ht, icr, itre, prob), sets imc=2, and computes crown width via cwcalc, so
   the example-tree machinery is properly populated and extree will not fault. Harness: read a keyword file
   that sets up the stand with NO DATABASE block (so no STANDSQL/TREESQL), run to stop point 7, call
   fvsAddTrees with the FIA tree list, then step the cycles with the stop-point-5 override. The one unknown
   to solve is stand-level setup without a database (site index, plot design, sampling weights) via
   keywords or the species/stand attribute API; that is the focused next implementation step.
2. Alternatively, debug the in-process DBS SQL generation: capture the exact statement the in-process DBS
   extension builds (the "unrecognized token: '" suggests a stray or doubled quote in the generated SQL),
   and fix the quoting in the DBS extension or the keyword so the SQLite read succeeds in process.
3. Cross-check against rFVS, which drives the same fvs API stop points and DBS reads from R; if rFVS loads
   the same SQLite stand cleanly, the gap is in how fvs2py sets up the database connection, not the engine.

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

## fvsAddTrees milestone (2026-06-18, evening) - the segfault is gone

Added an add_trees() method to fvs2py wrapping fvsAddTrees, and tested the database-free load path: a
keyword file with NO DATABASE block, run to stop point 7 (zero trees), then add_trees() to load the stand
in memory, then project. Result: add_trees loads the stand (50 synthetic trees, ntrees goes 0 -> 50) and
the engine projects WITHOUT the extree.f90 segfault. The documented Route A blocker is resolved: bypassing
the in-process DBS read via fvsAddTrees populates the tree arrays so extree no longer indexes unset
example-tree slots.

Remaining issue, downstream and more tractable: during PROCESS the reporting routine prtexm.f90 line 69
(print example trees) hits a Fortran I/O error, I/O past end of record on an unformatted scratch file
(unit 8), together with an OPEN FAILED FOR 17 (output database ref -1). The in-process API path has not
opened the main output and scratch files that the standalone executable opens, so the example-tree print
fails. Next: open or suppress the FVS output/scratch units in the in-process path (set the output file via
the API, or suppress the example-tree and main-output writes since results are read through the summary
API), then read the projected summary. This is output-unit plumbing, not a growth or memory fault.

Net: the in-memory tree load works and the segfault is eliminated. One output-init step remains before a
clean in-process projection, after which the fvs-conus DG/HD predictions wire into the stop-point-5 loop
for the true in-engine arms C and D.
