#!/bin/bash
# =============================================================================
# build_fvs_executables.sh
#
# Build FVS standalone executables for use with subprocess-based projection
# pipelines. Complements build_fvs_libraries.sh (which builds .so files).
#
# The key difference: this script links WITH main.o to produce a standalone
# binary that reads keyword files from stdin, rather than a shared library
# for ctypes/rFVS.
#
# Usage:
#   ./build_fvs_executables.sh [SOURCE_DIR] [OUTPUT_DIR] [VARIANTS...]
#
# Examples:
#   ./build_fvs_executables.sh . ./lib ne acd
#   ./build_fvs_executables.sh /path/to/fvs-modern /path/to/lib ne acd
#
# If no variants specified, builds ne and acd (the PERSEUS variants).
# =============================================================================

set -euo pipefail

SOURCE_DIR="${1:-.}"
OUTPUT_DIR="${2:-./lib}"
shift 2 2>/dev/null || true
VARIANTS=("$@")

# Default: PERSEUS variants
if [ ${#VARIANTS[@]} -eq 0 ]; then
    VARIANTS=(ne acd)
fi

# Compiler settings (match build_fvs_libraries.sh).
# gfortran is the only supported Fortran compiler. See the header comment in
# build_fvs_libraries.sh for the !DEC$ ATTRIBUTES ALIAS rationale.
if [ "${FC:-gfortran}" != "gfortran" ]; then
    echo "WARNING: FC=$FC is not supported. Reverting to gfortran." >&2
fi
FC="gfortran"
CC="${CC:-gcc}"
CXX="${CXX:-g++}"
# Note: no -fPIC needed for executables, but it doesn't hurt
FFLAGS="-g -cpp -DCMPgcc -std=legacy -w -fallow-argument-mismatch"
CFLAGS="-DANSI -DCMPgcc -DSQLITE_THREADSAFE=0 -DSQLITE_OMIT_LOAD_EXTENSION -w"
CXXFLAGS="-DANSI -DCMPgcc -w"

BUILD_DIR=$(mktemp -d /tmp/fvs-exe-build-XXXXXX)
trap "rm -rf $BUILD_DIR" EXIT

mkdir -p "$OUTPUT_DIR"

# Resolve source directory to converted source
if [ -d "$SOURCE_DIR/src-converted" ]; then
    SRC_ROOT=$(realpath "$SOURCE_DIR/src-converted")
elif [ -d "$SOURCE_DIR/base" ]; then
    SRC_ROOT=$(realpath "$SOURCE_DIR")
else
    echo "ERROR: Cannot find FVS source in $SOURCE_DIR"
    echo "Expected src-converted/ or base/ subdirectory."
    exit 1
fi

echo "================================================================"
echo "FVS Executable Builder"
echo "================================================================"
echo "Source:   $SRC_ROOT"
echo "Output:   $OUTPUT_DIR"
echo "Variants: ${VARIANTS[*]}"
echo "Compiler: $($FC --version | head -1)"
echo "Build:    $BUILD_DIR"
echo "================================================================"
echo ""

compile_file() {
    local src="$1"
    local obj="$2"
    local incdirs="$3"
    local ext="${src##*.}"

    case "$ext" in
        f|F|f90|F90|for)
            $FC $FFLAGS $incdirs -c "$src" -o "$obj" 2>/dev/null
            ;;
        c)
            $CC $CFLAGS $incdirs -c "$src" -o "$obj" 2>/dev/null
            ;;
        cpp|cxx|C)
            $CXX $CXXFLAGS $incdirs -c "$src" -o "$obj" 2>/dev/null
            ;;
    esac
}

BUILT=0
FAILED=0

for var in "${VARIANTS[@]}"; do
    SRCLIST="$SRC_ROOT/bin/FVS${var}_sourceList.txt"
    if [ ! -f "$SRCLIST" ]; then
        echo "[$var] SKIP: No source list at $SRCLIST"
        FAILED=$((FAILED + 1))
        continue
    fi

    echo -n "[$var] Building FVS${var} executable ... "
    VARDIR="$BUILD_DIR/$var"
    mkdir -p "$VARDIR"

    # 2026-05-16 holoros fork: shield variant-specific PRGPRM.f90 / ESPARM.f90
    # from being shadowed by base/ copies during compilation.
    #
    # gfortran resolves `INCLUDE 'PRGPRM.f90'` first against the current source
    # file's directory, THEN through -I dirs. Files in base/ (e.g. intree.f90,
    # keywds.f90) thus pick up base/PRGPRM.f90 (MAXSP=23) instead of the
    # variant's correct $var/common/PRGPRM.f90 (MAXSP=108 for ACD). The
    # resulting MAXSP discrepancy across compilation units corrupts the
    # /CONCHR/ common block layout: IUSED(MAXSP) sized at 4*23 bytes by base
    # files and 4*108 bytes by variant files puts TREFMT at different offsets
    # and the tree-format keyword handler writes through one offset while
    # intree.f90 reads through the other, garbling READ(RECORD,TREFMT) with
    # NAMGRP/PTGNAME data and aborting at intree.f90:188 with
    # "Missing initial left parenthesis in format".
    #
    # Move the base shadows aside for the duration of this variant build,
    # so all INCLUDE 'PRGPRM.f90' resolve via -I$var/common consistently.
    # Restore on exit so other variant builds and rFVS/fvsOL paths see the
    # original tree.
    SHADOWED_FILES=()
    for shadow_name in PRGPRM.f90 ESPARM.f90; do
        shadow_path="$SRC_ROOT/base/$shadow_name"
        variant_path="$SRC_ROOT/$var/common/$shadow_name"
        if [ -f "$shadow_path" ] && [ -f "$variant_path" ]; then
            mv "$shadow_path" "${shadow_path}.SHADOW_BUILD_${var}"
            SHADOWED_FILES+=("${shadow_path}")
        fi
    done
    restore_shadows() {
        for sp in "${SHADOWED_FILES[@]}"; do
            [ -f "${sp}.SHADOW_BUILD_${var}" ] && mv "${sp}.SHADOW_BUILD_${var}" "$sp"
        done
    }
    trap restore_shadows EXIT

    OBJECTS=()
    COMPILE_ERRORS=0

    while IFS= read -r srcfile || [[ -n "$srcfile" ]]; do
        srcfile=$(echo "$srcfile" | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        [ -z "$srcfile" ] && continue
        [[ "$srcfile" == \#* ]] && continue

        # Resolve relative paths from bin/ directory (source lists use ../
        # relative to their own location in bin/)
        if [[ ! "$srcfile" = /* ]]; then
            srcfile="$(cd "$SRC_ROOT/bin" && realpath -m "$srcfile")"
        fi

        if [ ! -f "$srcfile" ]; then
            # Source lists reference .f90 for files not yet converted from
            # fixed-form. Fall back to .f or .for.
            stem="${srcfile%.*}"
            if [ -f "${stem}.f" ]; then
                srcfile="${stem}.f"
            elif [ -f "${stem}.for" ]; then
                srcfile="${stem}.for"
            else
                COMPILE_ERRORS=$((COMPILE_ERRORS + 1))
                continue
            fi
        fi

        # Use unique object names to avoid collisions from different directories
        objname=$(echo "$srcfile" | sed "s|$SRC_ROOT/||" | tr '/' '_' | sed 's/\.[^.]*$/.o/')
        objpath="$VARDIR/$objname"

        # Include directories: variant common, base common, dbs, dbsqlite, fire, etc.
        # Note: base/ must come early so properly-converted .f90 includes are found
        INCDIRS="-I$SRC_ROOT/$var/common -I$SRC_ROOT/common -I$SRC_ROOT/base/common"
        INCDIRS="$INCDIRS -I$SRC_ROOT/$var -I$SRC_ROOT/base"
        INCDIRS="$INCDIRS -I$SRC_ROOT/dbs -I$SRC_ROOT/dbsqlite"
        [ -d "$SRC_ROOT/fire/$var" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/fire/$var"
        [ -d "$SRC_ROOT/fire/common" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/fire/common"
        [ -d "$SRC_ROOT/fire/base" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/fire/base"
        [ -d "$SRC_ROOT/fire/base/common" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/fire/base/common"
        [ -d "$SRC_ROOT/fire/$var/common" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/fire/$var/common"
        [ -d "$SRC_ROOT/estb" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/estb"
        [ -d "$SRC_ROOT/volume" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/volume"
        [ -d "$SRC_ROOT/volume/NVEL" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/volume/NVEL"
        [ -d "$SRC_ROOT/fiavbc" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/fiavbc"
        [ -d "$SRC_ROOT/vbase" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/vbase"
        [ -d "$SRC_ROOT/vvolume" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/vvolume"
        [ -d "$SRC_ROOT/vie" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/vie"
        [ -d "$SRC_ROOT/strp" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/strp"
        [ -d "$SRC_ROOT/clim/base" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/clim/base"
        [ -d "$SRC_ROOT/clim/$var" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/clim/$var"
        [ -d "$SRC_ROOT/vdbs" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/vdbs"
        [ -d "$SRC_ROOT/vdbsqlite" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/vdbsqlite"
        [ -d "$SRC_ROOT/econ" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/econ"
        [ -d "$SRC_ROOT/wpbr" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/wpbr"
        [ -d "$SRC_ROOT/wsbwe" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/wsbwe"
        [ -d "$SRC_ROOT/pg" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/pg"
        [ -d "$SRC_ROOT/canada/$var" ] && INCDIRS="$INCDIRS -I$SRC_ROOT/canada/$var"

        # Skip duplicate source list entries (same object name)
        if [ -f "$objpath" ]; then
            continue
        fi

        if compile_file "$srcfile" "$objpath" "$INCDIRS" 2>"$VARDIR/compile_${objname%.o}.err" && [ -f "$objpath" ]; then
            OBJECTS+=("$objpath")
        else
            # If .f90 failed to compile, try the .f fixed-form version
            # (common when .f90 has unconverted INCLUDE references)
            stem="${srcfile%.*}"
            fallback=""
            if [[ "$srcfile" == *.f90 ]] && [ -f "${stem}.f" ]; then
                fallback="${stem}.f"
            elif [[ "$srcfile" == *.f90 ]] && [ -f "${stem}.for" ]; then
                fallback="${stem}.for"
            fi
            if [ -n "$fallback" ]; then
                fb_objname=$(echo "$fallback" | sed "s|$SRC_ROOT/||" | tr '/' '_' | sed 's/\.[^.]*$/.o/')
                fb_objpath="$VARDIR/$fb_objname"
                if compile_file "$fallback" "$fb_objpath" "$INCDIRS" 2>/dev/null && [ -f "$fb_objpath" ]; then
                    OBJECTS+=("$fb_objpath")
                else
                    COMPILE_ERRORS=$((COMPILE_ERRORS + 1))
                fi
            else
                COMPILE_ERRORS=$((COMPILE_ERRORS + 1))
            fi
        fi

    done < "$SRCLIST"

    # Compile additional source directories that the source list references
    # only as .inc includes, not as compilable units. These provide symbols
    # like volinit2_, volinitnvb_, fmatv_, fmppget_, etc.
    EXTRA_SRC_DIRS=(
        "$SRC_ROOT/volume/NVEL"
    )
    EXTRA_INCDIRS="-I$SRC_ROOT/common -I$SRC_ROOT/$var/common -I$SRC_ROOT/base"
    EXTRA_INCDIRS="$EXTRA_INCDIRS -I$SRC_ROOT/fire/base -I$SRC_ROOT/fire/common -I$SRC_ROOT/fire/$var"
    EXTRA_INCDIRS="$EXTRA_INCDIRS -I$SRC_ROOT/volume/NVEL -I$SRC_ROOT/$var"

    for extra_dir in "${EXTRA_SRC_DIRS[@]}"; do
        [ -d "$extra_dir" ] || continue

        # Phase 1: Compile MODULE files first so .mod files are available
        # for subsequent compilations. Put .mod files in VARDIR.
        MODULE_INCDIRS="$EXTRA_INCDIRS -J$VARDIR -I$VARDIR"
        for mod_src in "$extra_dir"/clkcoef_mod.f "$extra_dir"/debug_mod.f "$extra_dir"/volinput_mod.f; do
            [ -f "$mod_src" ] || continue
            mod_base=$(basename "$mod_src")
            mod_obj="$VARDIR/extra_${mod_base%.*}.o"
            # Skip if already compiled from source list
            srclist_objname=$(echo "$mod_src" | sed "s|$SRC_ROOT/||" | tr '/' '_' | sed 's/\.[^.]*$/.o/')
            [ -f "$VARDIR/$srclist_objname" ] && continue
            if compile_file "$mod_src" "$mod_obj" "$MODULE_INCDIRS" 2>/dev/null && [ -f "$mod_obj" ]; then
                OBJECTS+=("$mod_obj")
            fi
        done

        # Phase 2: Compile all other source files.
        # Skip .f files when a .f90 free-form conversion exists (avoids double-compile).
        # Skip module files already compiled in Phase 1.
        PHASE2_INCDIRS="$EXTRA_INCDIRS -I$VARDIR"
        for extra_src in "$extra_dir"/*.f90 "$extra_dir"/*.f "$extra_dir"/*.for "$extra_dir"/*.c; do
            [ -f "$extra_src" ] || continue
            extra_base=$(basename "$extra_src")
            extra_stem="${extra_base%.*}"
            # Skip test files
            [[ "$extra_base" == *test* ]] && continue
            # Skip module files (already compiled)
            [[ "$extra_stem" == *_mod ]] && continue
            # Skip .f files when a .f90 conversion exists
            if [[ "$extra_base" == *.f ]] && [ -f "$extra_dir/${extra_stem}.f90" ]; then
                continue
            fi
            extra_obj="$VARDIR/extra_${extra_stem}.o"
            # Skip if already compiled (e.g., .f90 and .f both exist)
            [ -f "$extra_obj" ] && continue
            # Skip if already compiled from source list. Build the
            # expected object name the same way the source-list loop does
            # and check for that exact file.
            srclist_objname=$(echo "$extra_src" | sed "s|$SRC_ROOT/||" | tr '/' '_' | sed 's/\.[^.]*$/.o/')
            if [ -f "$VARDIR/$srclist_objname" ]; then
                continue
            fi
            if compile_file "$extra_src" "$extra_obj" "$PHASE2_INCDIRS" 2>/dev/null && [ -f "$extra_obj" ]; then
                OBJECTS+=("$extra_obj")
            fi
        done
    done

    # Compile stubs for subroutines that cannot compile on GCC
    # (Windows DLL imports, missing region-specific volume routines, etc.)
    # base/fvs_stubs.f90 is unconditional (stubs Windows DLL imports).
    for stub_src in "$SRC_ROOT/base/fvs_stubs.f90"; do
        [ -f "$stub_src" ] || continue
        stub_obj="$VARDIR/$(basename ${stub_src%.*}).o"
        if compile_file "$stub_src" "$stub_obj" "$INCDIRS" 2>/dev/null && [ -f "$stub_obj" ]; then
            OBJECTS+=("$stub_obj")
        fi
    done

    # Conditional stubs: only compile if the real implementation wasn't
    # successfully compiled from the source list. This handles cases like
    # fire/*/fmcfmd.f90 which uses LOC() that GCC rejects for some variants.
    HAS_FMCFMD=0
    for obj in "${OBJECTS[@]}"; do
        bobj=$(basename "$obj" .o)
        # Match fmcfmd but NOT fmcfmd2, fmcfmd3, fmcfmd_stub etc.
        case "$bobj" in *fmcfmd) HAS_FMCFMD=1; break;; esac
    done
    if [ "$HAS_FMCFMD" -eq 0 ] && [ -f "$SRC_ROOT/fire/vbase/fmcfmd_stub.f90" ]; then
        stub_obj="$VARDIR/fmcfmd_stub.o"
        if compile_file "$SRC_ROOT/fire/vbase/fmcfmd_stub.f90" "$stub_obj" "$INCDIRS" 2>/dev/null && [ -f "$stub_obj" ]; then
            OBJECTS+=("$stub_obj")
        fi
    fi

    # Conditional econ stubs (2026-05-16): same pattern as fmcfmd. The real
    # econ implementation (econ/ecin.f90, ecstatus.f90, eccalc.f90, etc.)
    # compiles cleanly under GCC for the ACD/NE/CS/LS/SN variants, so
    # econ_stubs.f90 would produce duplicate-symbol link errors when both
    # are present (eckey_ from ecin.f90:656 vs econ_stubs.f90:43;
    # getispretendactive_ from ecstatus.f90:80 vs econ_stubs.f90:48).
    # Only add the stubs file when the real econ objects are missing
    # (e.g., minimal-variant builds that exclude the econ tree).
    HAS_ECIN=0
    HAS_ECSTATUS=0
    for obj in "${OBJECTS[@]}"; do
        bobj=$(basename "$obj" .o)
        case "$bobj" in
            ecin|econ_ecin)         HAS_ECIN=1 ;;
            ecstatus|econ_ecstatus) HAS_ECSTATUS=1 ;;
        esac
    done
    if [ "$HAS_ECIN" -eq 0 ] && [ "$HAS_ECSTATUS" -eq 0 ] && [ -f "$SRC_ROOT/econ/econ_stubs.f90" ]; then
        stub_obj="$VARDIR/econ_stubs.o"
        if compile_file "$SRC_ROOT/econ/econ_stubs.f90" "$stub_obj" "$INCDIRS" 2>/dev/null && [ -f "$stub_obj" ]; then
            OBJECTS+=("$stub_obj")
        fi
    else
        # When econ_stubs.f90 is skipped (because real econ provides ecin /
        # ecstatus), we still need VARVER. It is called from variant code
        # (e.g. acd/dgdriv.f90:295, acd/morts.f90:139) but no upstream
        # variant actually defines it. Provide a tiny stub.
        if [ -f "$SRC_ROOT/base/varver_stub.f90" ]; then
            stub_obj="$VARDIR/varver_stub.o"
            if compile_file "$SRC_ROOT/base/varver_stub.f90" "$stub_obj" "$INCDIRS" 2>/dev/null && [ -f "$stub_obj" ]; then
                OBJECTS+=("$stub_obj")
            fi
        fi
    fi

    # Link executable (INCLUDING main.o, unlike the .so build)
    if [ ${#OBJECTS[@]} -gt 0 ]; then
        EXE_PATH="$OUTPUT_DIR/FVS${var}"
        # Try direct linking first
        if $FC -o "$EXE_PATH" "${OBJECTS[@]}" -lm 2>"$VARDIR/link.err"; then
            chmod +x "$EXE_PATH"
            NOBJ=${#OBJECTS[@]}
            SIZE=$(ls -lh "$EXE_PATH" | awk '{print $5}')
            echo "DONE ($NOBJ objects, $COMPILE_ERRORS skipped, $SIZE)"
            BUILT=$((BUILT + 1))
        else
            # Fallback: link against existing .so for missing symbols
            SO_PATH="$OUTPUT_DIR/FVS${var}.so"
            if [ -f "$SO_PATH" ]; then
                echo -n "direct link failed, trying with .so fallback ... "
                SO_DIR=$(dirname "$SO_PATH")
                SO_NAME="FVS${var}"
                if $FC -o "$EXE_PATH" "${OBJECTS[@]}" -L"$SO_DIR" -l:"${SO_NAME}.so" -Wl,-rpath,"$SO_DIR" -lm 2>"$VARDIR/link2.err"; then
                    chmod +x "$EXE_PATH"
                    NOBJ=${#OBJECTS[@]}
                    SIZE=$(ls -lh "$EXE_PATH" | awk '{print $5}')
                    echo "DONE via .so ($NOBJ objects, $COMPILE_ERRORS skipped, $SIZE)"
                    BUILT=$((BUILT + 1))
                else
                    echo "LINK FAILED (both direct and .so fallback)"
                    tail -5 "$VARDIR/link2.err" 2>/dev/null
                    FAILED=$((FAILED + 1))
                fi
            else
                echo "LINK FAILED (see $VARDIR/link.err)"
                tail -5 "$VARDIR/link.err" 2>/dev/null
                FAILED=$((FAILED + 1))
            fi
        fi
    else
        echo "NO OBJECTS"
        FAILED=$((FAILED + 1))
    fi

    # Restore base/ shadow files at end of this variant's iteration
    restore_shadows
    trap - EXIT
done

echo ""
echo "================================================================"
echo "Build Summary"
echo "================================================================"
echo "  Built successfully: $BUILT"
echo "  Failed:             $FAILED"
echo "  Output directory:   $OUTPUT_DIR"
echo ""
if [ $BUILT -gt 0 ]; then
    echo "Executables:"
    for var in "${VARIANTS[@]}"; do
        exe="$OUTPUT_DIR/FVS${var}"
        if [ -f "$exe" ]; then
            ls -lh "$exe"
        fi
    done
fi
echo "================================================================"
