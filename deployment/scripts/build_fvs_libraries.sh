#!/bin/bash
# =============================================================================
# build_fvs_libraries.sh
#
# Build FVS Fortran shared libraries (.so) for use with rFVS/fvsOL on Linux.
# Targets the ORIGINAL fixed-form source (ForestVegetationSimulator-main).
#
# Usage:
#   ./build_fvs_libraries.sh [--verbose] SOURCE_DIR OUTPUT_DIR [VARIANTS...]
#
# Examples:
#   ./build_fvs_libraries.sh . ./lib                    # all variants
#   ./build_fvs_libraries.sh . ./lib ne ak ie           # subset
#   ./build_fvs_libraries.sh --verbose . ./lib ne       # show compile errors
#
# If no variants specified, builds all US + Canadian variants.
# =============================================================================

set -euo pipefail

# Parse flags
VERBOSE=0
ARGS=()
for arg in "$@"; do
    case "$arg" in
        --verbose|-v) VERBOSE=1 ;;
        *) ARGS+=("$arg") ;;
    esac
done
set -- "${ARGS[@]+"${ARGS[@]}"}"

SOURCE_DIR="${1:?Usage: $0 [--verbose] SOURCE_DIR OUTPUT_DIR [VARIANTS...]}"
OUTPUT_DIR="${2:?Usage: $0 [--verbose] SOURCE_DIR OUTPUT_DIR [VARIANTS...]}"
shift 2
VARIANTS=("$@")

# Default: all US variants
if [ ${#VARIANTS[@]} -eq 0 ]; then
    VARIANTS=(ak acd bm ca ci cr cs ec em ie kt ls nc ne oc op pn sn so tt ut wc ws bc on)
fi

# Compiler settings
FC="${FC:-gfortran}"
CC="${CC:-gcc}"
CXX="${CXX:-g++}"
FFLAGS="-fPIC -g -cpp -DCMPgcc -std=legacy -w -fallow-argument-mismatch"
CFLAGS="-fPIC -DANSI -DCMPgcc -DSQLITE_THREADSAFE=0 -DSQLITE_OMIT_LOAD_EXTENSION -w"
CXXFLAGS="-fPIC -DANSI -DCMPgcc -w"

BUILD_DIR=$(mktemp -d /tmp/fvs-build-XXXXXX)
trap "rm -rf $BUILD_DIR" EXIT

mkdir -p "$OUTPUT_DIR"

echo "================================================================"
echo "FVS Shared Library Builder"
echo "================================================================"
echo "Source:   $SOURCE_DIR"
echo "Output:   $OUTPUT_DIR"
echo "Variants: ${VARIANTS[*]}"
echo "Compiler: $($FC --version | head -1)"
echo "Build:    $BUILD_DIR"
echo "================================================================"
echo ""

# Resolve source directory: accept project root (with src-converted/) or
# direct path to the source tree (with bin/ and base/).
if [ -d "$SOURCE_DIR/src-converted/bin" ] && [ -d "$SOURCE_DIR/src-converted/base" ]; then
    SOURCE_DIR="$SOURCE_DIR/src-converted"
elif [ ! -d "$SOURCE_DIR/bin" ] || [ ! -d "$SOURCE_DIR/base" ]; then
    echo "ERROR: $SOURCE_DIR does not look like a FVS source tree."
    echo "Expected bin/ and base/ subdirectories (or src-converted/ parent)."
    exit 1
fi

# Function to compile a single source file
compile_file() {
    local src="$1"
    local obj="$2"
    local incdirs="$3"
    local ext="${src##*.}"

    case "$ext" in
        f|F|f90|F90)
            $FC $FFLAGS $incdirs -c "$src" -o "$obj" 2>/dev/null
            ;;
        for)
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

# Build each variant
BUILT=0
FAILED=0

for var in "${VARIANTS[@]}"; do
    SRCLIST="$SOURCE_DIR/bin/FVS${var}_sourceList.txt"
    if [ ! -f "$SRCLIST" ]; then
        echo "[$var] SKIP: No source list found at $SRCLIST"
        continue
    fi

    echo -n "[$var] Building FVS${var}.so ... "
    VARDIR="$BUILD_DIR/$var"
    mkdir -p "$VARDIR"

    # Parse source list and collect include directories
    OBJECTS=()
    COMPILE_ERRORS=0

    while IFS= read -r line || [ -n "$line" ]; do
        # Skip blank lines and comments
        line=$(echo "$line" | sed 's/^[[:space:]]*//' | sed 's/[[:space:]]*$//')
        [ -z "$line" ] && continue
        [[ "$line" == "#"* ]] && continue

        # Resolve path relative to bin/
        srcfile="$SOURCE_DIR/bin/$line"
        if [ ! -f "$srcfile" ]; then
            # Source lists often reference .f90 for files that haven't been
            # converted from fixed-form yet. Fall back to .f or .for.
            stem="${srcfile%.*}"
            if [ -f "${stem}.f" ]; then
                srcfile="${stem}.f"
            elif [ -f "${stem}.for" ]; then
                srcfile="${stem}.for"
            else
                continue
            fi
        fi

        # Build include path from source directory
        srcdir=$(dirname "$srcfile")
        fname=$(basename "$srcfile")
        objname="${fname%.*}.o"
        objpath="$VARDIR/$objname"

        # Handle duplicate basenames by adding directory prefix
        if [ -f "$objpath" ]; then
            dirprefix=$(basename "$srcdir")
            objpath="$VARDIR/${dirprefix}_${objname}"
        fi

        # Collect include directories. Variant-specific paths come BEFORE
        # base/ so that PRGPRM.f90 and other shared header names resolve
        # to the variant's own copy (with the variant-specific MAXSP and
        # other parameters). Without this, files like vwc/morts.f90 that
        # are shared between PN and WC and depend on MAXSP=39 would
        # silently pick up base/PRGPRM.f90's MAXSP=23 and fail to
        # compile their 39-element DATA arrays.
        INCDIRS="-I$SOURCE_DIR/$var"
        [ -d "$SOURCE_DIR/$var/common" ] && INCDIRS="$INCDIRS -I$SOURCE_DIR/$var/common"
        INCDIRS="$INCDIRS -I$SOURCE_DIR/base"
        [ -d "$SOURCE_DIR/base/common" ] && INCDIRS="$INCDIRS -I$SOURCE_DIR/base/common"
        [ -d "$SOURCE_DIR/vbase" ] && INCDIRS="$INCDIRS -I$SOURCE_DIR/vbase"
        [ -d "$SOURCE_DIR/common" ] && INCDIRS="$INCDIRS -I$SOURCE_DIR/common"
        [ -d "$SOURCE_DIR/fire/base/common" ] && INCDIRS="$INCDIRS -I$SOURCE_DIR/fire/base/common"
        [ -d "$SOURCE_DIR/fire/$var/common" ] && INCDIRS="$INCDIRS -I$SOURCE_DIR/fire/$var/common"
        # Canadian variants live under canada/ subdirectory
        [ -d "$SOURCE_DIR/canada/$var" ] && INCDIRS="$INCDIRS -I$SOURCE_DIR/canada/$var"
        [ -d "$SOURCE_DIR/canada/fire" ] && INCDIRS="$INCDIRS -I$SOURCE_DIR/canada/fire"
        [ -d "$SOURCE_DIR/canada/newmist" ] && INCDIRS="$INCDIRS -I$SOURCE_DIR/canada/newmist"
        [ -d "$srcdir" ] && INCDIRS="$INCDIRS -I$srcdir"
        INCDIRS="$INCDIRS -I$SOURCE_DIR/dbs -I$SOURCE_DIR/dbsqlite"
        # NVEL volume library includes (.inc files like wdbkwtdata.inc)
        [ -d "$SOURCE_DIR/volume/NVEL" ] && INCDIRS="$INCDIRS -I$SOURCE_DIR/volume/NVEL"
        [ -d "$SOURCE_DIR/volume" ] && INCDIRS="$INCDIRS -I$SOURCE_DIR/volume"
        # Cross-variant source files (e.g., ie/vols.f90 used by acd) need
        # the originating variant's common dir for ESPARM.f90 etc.
        [ -d "$SOURCE_DIR/ie/common" ] && INCDIRS="$INCDIRS -I$SOURCE_DIR/ie/common"
        [ -d "$SOURCE_DIR/ls/common" ] && INCDIRS="$INCDIRS -I$SOURCE_DIR/ls/common"

        if compile_file "$srcfile" "$objpath" "$INCDIRS" 2>"$VARDIR/err_$(basename "$objpath" .o).log" && [ -f "$objpath" ]; then
            OBJECTS+=("$objpath")
        else
            # If .f90 failed, try .f fixed-form fallback (common when .f90
            # has unconverted INCLUDE references but .f compiles cleanly)
            stem="${srcfile%.*}"
            fallback=""
            if [[ "$srcfile" == *.f90 ]] && [ -f "${stem}.f" ]; then
                fallback="${stem}.f"
            elif [[ "$srcfile" == *.f90 ]] && [ -f "${stem}.for" ]; then
                fallback="${stem}.for"
            fi
            if [ -n "$fallback" ]; then
                fb_objname="${fallback##*/}"
                fb_objpath="$VARDIR/fb_${fb_objname%.*}.o"
                if compile_file "$fallback" "$fb_objpath" "$INCDIRS" 2>/dev/null && [ -f "$fb_objpath" ]; then
                    OBJECTS+=("$fb_objpath")
                else
                    COMPILE_ERRORS=$((COMPILE_ERRORS + 1))
                    if [ "$VERBOSE" -eq 1 ]; then
                        echo ""
                        echo "  COMPILE FAIL: $srcfile (and fallback $fallback)"
                        tail -3 "$VARDIR/err_$(basename "$objpath" .o).log" 2>/dev/null | sed 's/^/    /'
                    fi
                fi
            else
                COMPILE_ERRORS=$((COMPILE_ERRORS + 1))
                if [ "$VERBOSE" -eq 1 ]; then
                    echo ""
                    echo "  COMPILE FAIL: $srcfile"
                    tail -3 "$VARDIR/err_$(basename "$objpath" .o).log" 2>/dev/null | sed 's/^/    /'
                fi
            fi
        fi

    done < "$SRCLIST"

    # Filter out main.o (the executable entry point; not needed for shared library)
    SHLIB_OBJECTS=()
    for obj in "${OBJECTS[@]}"; do
        bname=$(basename "$obj")
        if [[ "$bname" != "main.o" ]]; then
            SHLIB_OBJECTS+=("$obj")
        fi
    done

    # Link shared library
    if [ ${#SHLIB_OBJECTS[@]} -gt 0 ]; then
        if $FC -shared -o "$OUTPUT_DIR/FVS${var}.so" "${SHLIB_OBJECTS[@]}" 2>"$VARDIR/link.err"; then
            NOBJ=${#SHLIB_OBJECTS[@]}
            SIZE=$(ls -lh "$OUTPUT_DIR/FVS${var}.so" | awk '{print $5}')
            echo "DONE ($NOBJ objects, $COMPILE_ERRORS skipped, $SIZE)"
            BUILT=$((BUILT + 1))
        else
            echo "LINK FAILED"
            if [ "$VERBOSE" -eq 1 ]; then
                tail -5 "$VARDIR/link.err" 2>/dev/null | sed 's/^/    /'
            fi
            FAILED=$((FAILED + 1))
        fi
    else
        echo "NO OBJECTS"
        FAILED=$((FAILED + 1))
    fi
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
    echo "Libraries:"
    ls -lh "$OUTPUT_DIR"/FVS*.so 2>/dev/null
fi
echo "================================================================"
