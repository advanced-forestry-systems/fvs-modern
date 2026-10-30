# Building the engine

The engine is built as one shared library per variant, `FVS<variant>.so`, by `deployment/scripts/build_fvs_libraries.sh`. This page documents the script as it parses its arguments and as it behaves, not as older notes describe it.

## Toolchain

gfortran is the only supported Fortran compiler. The script forces `FC=gfortran`; if the environment sets `FC` to anything else it prints a warning and reverts. The stated reason is that Intel ifort exports uppercase symbols for routines decorated with `!DEC$ ATTRIBUTES ALIAS` (for example in `src-converted/base/fvs.f90`), which do not match the lowercase, underscore-suffixed names that Fortran callers and ctypes consumers expect. C sources are compiled with `${CC:-gcc}` and C++ sources with `${CXX:-g++}`.

`CONTRIBUTING.md` asks for gfortran 12 or later. On OSC Cardinal the script header notes that a session loads `gcc/12.3.0` first. Platform packages named in the script header are `gcc-gfortran` (Fedora/RHEL), `gfortran` (Ubuntu/Debian) and Homebrew `gcc` (macOS).

The compiler flags are fixed in the script:

```bash
FFLAGS="-fPIC -g -cpp -DCMPgcc -std=legacy -w -fallow-argument-mismatch"
CFLAGS="-fPIC -DANSI -DCMPgcc -DSQLITE_THREADSAFE=0 -DSQLITE_OMIT_LOAD_EXTENSION -w"
CXXFLAGS="-fPIC -DANSI -DCMPgcc -w"
```

`-w` suppresses all warnings, so a clean console does not mean a clean compile. Use `--verbose` and the per-file error logs described below.

## Usage

```bash
./build_fvs_libraries.sh [--verbose] SOURCE_DIR OUTPUT_DIR [VARIANTS...]
```

The script scans all arguments for `--verbose` or `-v` (anywhere on the line) and treats every other argument as positional. `SOURCE_DIR` and `OUTPUT_DIR` are required. Any further arguments are variant codes.

`SOURCE_DIR` may be the repository root or the source tree itself. If `SOURCE_DIR/src-converted/bin` and `SOURCE_DIR/src-converted/base` exist, the script descends into `src-converted`; otherwise it requires `bin/` and `base/` directly under `SOURCE_DIR` and exits with an error if they are missing. Both forms appear in the CI workflows:

```bash
bash deployment/scripts/build_fvs_libraries.sh src-converted lib ne
bash deployment/scripts/build_fvs_libraries.sh . lib ne ls
```

`--verbose` prints the failing file and the last three lines of its error log for each compile failure. It is a command-line flag only. The script initialises `VERBOSE=0` before parsing, so exporting a `VERBOSE` environment variable has no effect.

The script compiles into a temporary directory created with `mktemp -d /tmp/fvs-build-XXXXXX` and removes it on exit. Libraries are written to `OUTPUT_DIR`.

## Variants

With no variant arguments the script builds this list, in this order:

```text
ak acd bm ca ci cr cs ec em ie kt ls nc ne oc op pn sn so tt ut wc ws bc on
```

That is 25 variants, including the two Canadian variants `bc` and `on`. A variant is built only if `src-converted/bin/FVS<variant>_sourceList.txt` exists; otherwise the script prints `[<variant>] SKIP: No source list found` and moves on. The `bin/` directory also carries `FVSadk_sourceList.txt`; `adk` is not in the default list and must be named explicitly.

## What the script does per variant

1. **Stub libraries.** Before any variant, `build_stubs` compiles each `src-converted/stubs/*.f` into `OUTPUT_DIR/lib<name>.so` (currently `fvs_stubs.f` and `fvs_stubs_final.f`). On Linux the variant link adds these with `-l:` so they resolve only symbols the variant does not define, and sets an `$ORIGIN` rpath so the variant library finds them beside itself.
2. **PRGPRM shadow.** If the variant's `PRGPRM.f90` differs from `base/PRGPRM.f90`, the script copies `base/PRGPRM.f90` to `base/PRGPRM.f90.prgprm_bak`, overwrites `base/PRGPRM.f90` with the variant copy, and restores it at the end of the iteration and from the EXIT trap. The reason is that gfortran resolves a bare `INCLUDE 'PRGPRM.f90'` in a `base/` source against the source file's own directory first, regardless of `-I` order, so `base/` objects would otherwise be sized with `MAXSP=23` while variant objects use the variant value, corrupting `COMMON /CONCHR/`.
3. **Source list.** Each non-blank, non-`#` line of the source list is resolved relative to `bin/`. Five INCLUDE-only files under `rd/` (`RDADD`, `RDARRY`, `RDCOM`, `RDCRY`, `RDPARM`) are skipped. If a listed `.f90` does not exist, the script falls back to `.f` and then `.for`; if a `.f90` exists but fails to compile, it retries the `.f` or `.for` sibling.
4. **Include path.** Variant directories come before `base/` (`-I<var>`, `-I<var>/common`, then `-Ibase`, `-Ibase/common`, `-Ivbase`, `-Icommon`, fire, Canadian, `dbs`, `dbsqlite`, `volume/NVEL`, `volume`, `ie/common`, `ls/common`).
5. **Link.** `main.o` is excluded and the remaining objects are linked with `gfortran -shared`. Linux adds the stub libraries and rpath; macOS adds `-Wl,-undefined,dynamic_lookup`; MinGW/MSYS/Cygwin statically bundles the gfortran runtime and first generates pre-link stubs for undefined symbols.
6. **CONCHR assertion.** After linking, the script compares the size of the `conchr_` symbol in `intree.o` and `blkdat.o` using `nm --print-size`. A mismatch deletes the library, counts the variant as failed and prints `[CONCHR] ASSERTION FAILED`. A match prints `[CONCHR] OK`.
7. **Self-contained pass.** Remaining undefined lowercase Fortran symbols in the library are written as empty subroutines, compiled and relinked so the library also loads where the loader cannot defer unresolved symbols.

A successful variant prints a line of the form `DONE (<objects> objects, <skipped> skipped, <size>)`. The summary at the end reports `Built successfully` and `Failed` counts.

!!! warning "A DONE line can hide missing objects"
    A file that fails to compile (and has no working fixed-form fallback) is counted in the `skipped` figure and left out of the link. The library still links and the variant is still reported as built. Always read the skipped count, and rerun with `--verbose` or read `err_<object>.log` in the build directory when it is not what you expect for that variant.

## Build serially

The PRGPRM shadow edits `base/PRGPRM.f90` in place inside `SOURCE_DIR`. Two builds of different variants in the same source tree at the same time therefore see each other's `PRGPRM.f90`. Parallel variant builds in one tree are unsafe: a parallel build once dropped 7 objects from the NE library. Build variants serially, either in one invocation that lists several variants or in successive invocations. If you need parallelism, give each build its own copy of the source tree.

`deployment/scripts/build_fvs_executables.sh` has the same property. It moves `base/PRGPRM.f90` and `base/ESPARM.f90` aside to `*.SHADOW_BUILD_<var>` for the duration of a variant build, so it must not run concurrently with another build in the same tree either.

```bash
# One serial invocation for several variants
bash deployment/scripts/build_fvs_libraries.sh src-converted ./lib ne ak ie
```

## Standalone executables

`deployment/scripts/build_fvs_executables.sh [SOURCE_DIR] [OUTPUT_DIR] [VARIANTS...]` links `main.o` to produce standalone binaries. With no variants it builds `ne` and `acd`. The CI job that uses it first builds the shared libraries into the same output directory for fallback linking:

```bash
bash deployment/scripts/build_fvs_libraries.sh . lib ne ls
bash deployment/scripts/build_fvs_executables.sh . lib ne ls
```

## Verifying a build

Check three things: the build log, the exported API symbols, and that the library loads.

**Build log.** Each variant should print `[CONCHR] OK` and a `DONE` line whose skipped count you expect. A `[self-contained] stubbed N undefined symbol(s)` line is normal.

**Symbols.** The CI job checks for the API entry points with `nm` and fails if `fvssetcmdline_` or `fvssummary_` is missing. It also fails if `FVSne.so` is smaller than 1,000,000 bytes.

```bash
nm lib/FVSne.so | grep " [Tt] fvssetcmdline_"
```

**Load test with full symbol resolution.** `deployment/scripts/test_all_variants.sh` builds the requested variants and then loads every `FVS*.so`, `FVS*.dylib` or `FVS*.dll` in the output directory through ctypes with `RTLD_GLOBAL | RTLD_NOW`, so an undefined symbol fails the load instead of surfacing later. It exits non-zero if any library fails.

```bash
bash deployment/scripts/test_all_variants.sh src-converted lib ne
```

The core of its load step is:

```python
ctypes.CDLL(L, mode=getattr(ctypes,"RTLD_GLOBAL",0)|os.RTLD_NOW if hasattr(os,"RTLD_NOW") else 0)
```

Note that `test_all_variants.sh` invokes the library builder itself, so do not run it alongside another build in the same tree. The CI workflows load with `RTLD_LAZY` instead (see [Testing and CI](testing-and-ci.md)), which is a weaker check.

## Common failures

| Symptom | Cause in the script | Action |
|---|---|---|
| `WARNING: FC=... is not supported. Reverting to gfortran.` | `FC` set to a non-gfortran compiler | Informational; the build proceeds with gfortran. |
| `[<var>] SKIP: No source list found at ...` | No `bin/FVS<var>_sourceList.txt` | Check the variant code and `SOURCE_DIR`. |
| `ERROR: ... does not look like a FVS source tree.` | `SOURCE_DIR` has neither `src-converted/bin` plus `src-converted/base` nor `bin/` plus `base/` | Pass the repository root or `src-converted`. |
| High `skipped` count in the `DONE` line | Compile failures dropped from the link | Rerun with `--verbose`; check for a concurrent build in the same tree. |
| `[CONCHR] ASSERTION FAILED` | `intree.o` and `blkdat.o` disagree on `COMMON /CONCHR/` size, the zero-trees failure mode | The library is deleted. Check that `base/PRGPRM.f90` was not left shadowed by an interrupted or concurrent build. |
| `LINK FAILED` followed by the tail of `link.err` | Link error | Read the printed tail. |
| `NO OBJECTS` | Every file failed to compile | Rerun with `--verbose`. |
| `[stubs] Building ... FAIL (see /tmp/stub_<name>.err)` | Stub source failed to compile | Read the named error file; Linux variant links rely on these libraries for fallback resolution. |
| `base/PRGPRM.f90.prgprm_bak` present after a build | The build was killed before the restore ran | Restore `base/PRGPRM.f90` from the backup or from version control before building again. |

!!! note "Stale script comments"
    The script header still describes its target as the original fixed-form source tree and labels the default list as "all US variants". The code resolves `src-converted/` and the default list includes `bc` and `on`. The header text is stale; the behavior described on this page is taken from the code.
