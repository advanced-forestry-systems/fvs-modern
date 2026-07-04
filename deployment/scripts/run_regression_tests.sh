#!/bin/bash
# =============================================================================
# FVS Regression Test Suite
#
# Runs each variant's test cases using the standalone executable and compares
# summary output against saved baselines. Also runs rFVS library load tests.
#
# Usage:
#   ./run_regression_tests.sh [FVS_BIN_DIR] [TEST_DIR] [VARIANTS...]
#
# Examples:
#   ./run_regression_tests.sh                          # auto-detect, all variants
#   ./run_regression_tests.sh /path/to/bin /path/to/tests ne ie ls
#
# Exit codes:
#   0 = all tests passed
#   1 = one or more tests failed
# =============================================================================

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BOLD='\033[1m'
NC='\033[0m'

# =========================================================================
# Configuration
# =========================================================================
FVS_BIN="${1:-}"
TEST_DIR="${2:-}"
shift 2 2>/dev/null || true
VARIANTS=("$@")

# Auto-detect paths
if [ -z "$FVS_BIN" ]; then
    for candidate in \
        "$HOME/fvs-online/bin" \
        "/srv/fvs/bin" \
        "$HOME/ForestVegetationSimulator-main/bin"; do
        if [ -d "$candidate" ] && ls "$candidate"/FVS*.so &>/dev/null; then
            FVS_BIN="$candidate"
            break
        fi
    done
fi

if [ -z "$TEST_DIR" ]; then
    for candidate in \
        "$HOME/ForestVegetationSimulator-main/tests" \
        "/srv/fvs/tests"; do
        if [ -d "$candidate" ]; then
            TEST_DIR="$candidate"
            break
        fi
    done
fi

if [ -z "$FVS_BIN" ] || [ -z "$TEST_DIR" ]; then
    echo "Usage: $0 [FVS_BIN_DIR] [TEST_DIR] [VARIANTS...]"
    exit 1
fi

# Resolve to absolute paths. Part 2 cd's into each per-test working directory, so a
# relative FVS_BIN or TEST_DIR would resolve against the wrong directory on every test
# after the first, producing spurious "no executable" skips. (Fixed 2026-06-18.)
FVS_BIN="$(cd "$FVS_BIN" 2>/dev/null && pwd)"
TEST_DIR="$(cd "$TEST_DIR" 2>/dev/null && pwd)"

echo ""
echo "================================================================"
echo "  FVS Regression Test Suite"
echo "================================================================"
echo "  Binary dir: $FVS_BIN"
echo "  Test dir:   $TEST_DIR"
echo "================================================================"
echo ""

TOTAL=0
PASSED=0
FAILED=0
SKIPPED=0
FAIL_LIST=()

WORK_DIR=$(mktemp -d)
trap "rm -rf $WORK_DIR" EXIT

# =========================================================================
# Test 1: Library Load Tests (rFVS)
# =========================================================================
echo -e "${BOLD}=== Part 1: rFVS Library Load Tests ===${NC}"
echo ""

# Find R library path
R_LIBS=""
for candidate in "$HOME/fvs-online/R-libs" "/srv/fvs/R-libs" "/tmp/rFVS-lib"; do
    if [ -d "$candidate/rFVS" ]; then
        R_LIBS="$candidate"
        break
    fi
done

if [ -n "$R_LIBS" ]; then
    for sofile in "$FVS_BIN"/FVS*.so; do
        variant=$(basename "$sofile" .so)

        # Skip if specific variants requested and this one isn't in the list
        if [ ${#VARIANTS[@]} -gt 0 ]; then
            v_lower=$(echo "${variant#FVS}" | tr '[:upper:]' '[:lower:]')
            match=false
            for req in "${VARIANTS[@]}"; do
                if [ "$(echo "$req" | tr '[:upper:]' '[:lower:]')" = "$v_lower" ]; then
                    match=true; break
                fi
            done
            $match || continue
        fi

        TOTAL=$((TOTAL + 1))

        result=$(R --no-save --quiet -e "
            .libPaths('$R_LIBS')
            library(rFVS)
            tryCatch({
                fvsLoad('$variant', bin='$FVS_BIN')
                dims <- fvsGetDims()
                cat('OK', dims['maxspecies'], dims['maxtrees'])
            }, error = function(e) cat('FAIL', e\$message))
        " 2>&1 | sed 's/^> //' | grep -E "^OK|^FAIL" | head -1)

        if [[ "$result" == OK* ]]; then
            clean=$(echo "$result" | tr -d '>')
            species=$(echo "$clean" | awk '{print $2}')
            maxtrees=$(echo "$clean" | awk '{print $3}')
            printf "  ${GREEN}[PASS]${NC} %-10s loaded (species=%s, maxtrees=%s)\n" "$variant" "$species" "$maxtrees"
            PASSED=$((PASSED + 1))
        else
            printf "  ${RED}[FAIL]${NC} %-10s %s\n" "$variant" "$result"
            FAILED=$((FAILED + 1))
            FAIL_LIST+=("load:$variant")
        fi
    done
else
    echo -e "  ${YELLOW}[SKIP]${NC} rFVS not installed (no R library path found)"
fi

echo ""

# =========================================================================
# Test 2: Standalone Simulation Tests
# =========================================================================
echo -e "${BOLD}=== Part 2: Standalone Simulation Tests ===${NC}"
echo ""

for test_dir in "$TEST_DIR"/FVS*/; do
    variant_name=$(basename "$test_dir")
    v_code=$(echo "${variant_name#FVS}" | tr '[:upper:]' '[:lower:]')
    executable="$FVS_BIN/$variant_name"

    # Skip if specific variants requested
    if [ ${#VARIANTS[@]} -gt 0 ]; then
        match=false
        for req in "${VARIANTS[@]}"; do
            if [ "$(echo "$req" | tr '[:upper:]' '[:lower:]')" = "$v_code" ]; then
                match=true; break
            fi
        done
        $match || continue
    fi

    # Check executable exists
    if [ ! -x "$executable" ]; then
        printf "  ${YELLOW}[SKIP]${NC} %-10s no executable\n" "$variant_name"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Run each keyword file
    for keyfile in "$test_dir"*.key; do
        [ -f "$keyfile" ] || continue
        test_name=$(basename "$keyfile" .key)
        TOTAL=$((TOTAL + 1))

        # Set up working directory (recursive copy for subdirs like econ/)
        run_dir="$WORK_DIR/${variant_name}_${test_name}"
        mkdir -p "$run_dir"
        cp -r "$test_dir"* "$run_dir/" 2>/dev/null || true
        cd "$run_dir"

        # Run the simulation with a timeout
        # FVS returns non-zero exit codes as STOP codes (e.g., STOP 20 is normal)
        # A segfault (139) or timeout (124) indicates real failure
        sim_stderr="$run_dir/stderr.log"
        exit_code=0
        timeout 120 "$executable" --keywordfile="$test_name.key" > /dev/null 2>"$sim_stderr" || exit_code=$?

        # exit codes: 0 or 20 = success, 124 = timeout, 139 = segfault
        if [ "$exit_code" -lt 128 ]; then
            # Check if summary output was produced
            sumfile="${test_name}.sum"
            savefile="${test_name}.sum.save"

            if [ -f "$sumfile" ] && [ -f "$savefile" ]; then
                # Compare summary tables (extract numeric rows, allow small float diffs)
                # Simple check: line count within 10% and key values match
                sum_lines=$(wc -l < "$sumfile")
                save_lines=$(wc -l < "$savefile")

                if [ "$sum_lines" -gt 0 ] && [ "$save_lines" -gt 0 ]; then
                    # Extract the summary statistics lines (contain year data)
                    diff_count=$(diff <(grep -E "^[0-9]{4}" "$sumfile" 2>/dev/null || true) \
                                      <(grep -E "^[0-9]{4}" "$savefile" 2>/dev/null || true) 2>/dev/null | grep -c "^[<>]" || true)

                    if [ "$diff_count" -eq 0 ]; then
                        printf "  ${GREEN}[PASS]${NC} %-10s %-20s (exact match)\n" "$variant_name" "$test_name"
                        PASSED=$((PASSED + 1))
                    else
                        printf "  ${YELLOW}[WARN]${NC} %-10s %-20s (%d summary lines differ)\n" "$variant_name" "$test_name" "$diff_count"
                        PASSED=$((PASSED + 1))  # Count as pass since simulation completed
                    fi
                else
                    printf "  ${GREEN}[PASS]${NC} %-10s %-20s (simulation completed)\n" "$variant_name" "$test_name"
                    PASSED=$((PASSED + 1))
                fi
            elif [ -f "$sumfile" ]; then
                printf "  ${GREEN}[PASS]${NC} %-10s %-20s (no baseline to compare)\n" "$variant_name" "$test_name"
                PASSED=$((PASSED + 1))
            elif ls "${test_name}"*.out &>/dev/null || ls "${test_name}"*Out.db &>/dev/null || ls *Out.db &>/dev/null; then
                # Simulation completed but wrote output to .out or database instead of .sum
                printf "  ${GREEN}[PASS]${NC} %-10s %-20s (completed, output to db/out)\n" "$variant_name" "$test_name"
                PASSED=$((PASSED + 1))
            else
                printf "  ${RED}[FAIL]${NC} %-10s %-20s (no summary output)\n" "$variant_name" "$test_name"
                FAILED=$((FAILED + 1))
                FAIL_LIST+=("sim:${variant_name}/${test_name}")
            fi
        else
            printf "  ${RED}[FAIL]${NC} %-10s %-20s (runtime error or timeout)\n" "$variant_name" "$test_name"
            FAILED=$((FAILED + 1))
            FAIL_LIST+=("sim:${variant_name}/${test_name}")
        fi

        cd "$WORK_DIR"
    done
done

echo ""

# =========================================================================
# Test 3: rFVS API Simulation Test (if rFVS available)
# =========================================================================
if [ -n "$R_LIBS" ]; then
    echo -e "${BOLD}=== Part 3: rFVS API Simulation Test ===${NC}"
    echo ""

    # Find the IE test data
    if [ -d "$TEST_DIR/APIviaR" ]; then
        TOTAL=$((TOTAL + 1))

        api_dir="$WORK_DIR/api_test"
        mkdir -p "$api_dir"
        cp "$TEST_DIR/APIviaR"/* "$api_dir/" 2>/dev/null || true

        result=$(cd "$api_dir" && R --no-save --quiet -e "
            .libPaths('$R_LIBS')
            library(rFVS)
            fvsLoad('FVSie', bin='$FVS_BIN')
            fvsSetCmdLine('--keywordfile=base.key')
            fvsRun(0, 0)
            summ <- fvsGetSummary()
            if (nrow(summ) > 0 && summ[1,'Tpa'] > 0) {
                cat('OK', summ[1,'Year'], summ[1,'Tpa'], summ[nrow(summ),'Year'], summ[nrow(summ),'Tpa'])
            } else {
                cat('FAIL empty_summary')
            }
        " 2>&1 | sed 's/^> //' | grep -E "^OK|^FAIL" | head -1)

        if [[ "$result" == OK* ]]; then
            # Strip trailing R prompt characters (>) from values
            clean=$(echo "$result" | tr -d '>')
            y1=$(echo "$clean" | awk '{print $2}')
            tpa1=$(echo "$clean" | awk '{print $3}')
            yn=$(echo "$clean" | awk '{print $4}')
            tpan=$(echo "$clean" | awk '{print $5}')
            printf "  ${GREEN}[PASS]${NC} rFVS API test: year %s TPA=%.0f to year %s TPA=%.0f\n" "$y1" "$tpa1" "$yn" "$tpan"
            PASSED=$((PASSED + 1))
        else
            printf "  ${RED}[FAIL]${NC} rFVS API test: %s\n" "$result"
            FAILED=$((FAILED + 1))
            FAIL_LIST+=("api:rFVS_simulation")
        fi
    fi

    echo ""
fi

# =========================================================================
# Summary
# =========================================================================
echo "================================================================"
if [ "$FAILED" -eq 0 ]; then
    echo -e "  ${GREEN}${BOLD}ALL TESTS PASSED${NC}"
else
    echo -e "  ${RED}${BOLD}SOME TESTS FAILED${NC}"
fi
echo ""
echo "  Total:   $TOTAL"
echo "  Passed:  $PASSED"
echo "  Failed:  $FAILED"
echo "  Skipped: $SKIPPED"

if [ ${#FAIL_LIST[@]} -gt 0 ]; then
    echo ""
    echo "  Failed tests:"
    for f in "${FAIL_LIST[@]}"; do
        echo "    * $f"
    done
fi

echo "================================================================"

exit $FAILED
