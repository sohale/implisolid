#!/usr/bin/env bash
# Run all slang/native sanity tests in order.
# Stops on first failure unless --keep-going is passed.
#
# Usage:
#   bash scripts/run-all-tests.sh
#   bash scripts/run-all-tests.sh --keep-going

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KEEP_GOING=0
[[ "${1:-}" == "--keep-going" ]] && KEEP_GOING=1

PASS=0
FAIL=0

run_test() {
    local name="$1"; local script="$SCRIPT_DIR/$2"
    echo ""
    echo "════════════════════════════════════════"
    echo " $name"
    echo "════════════════════════════════════════"
    if bash "$script"; then
        echo "[SUITE PASS] $name"
        ((PASS++)) || true
    else
        echo "[SUITE FAIL] $name"
        ((FAIL++)) || true
        [[ $KEEP_GOING -eq 0 ]] && exit 1
    fi
}

run_test "IImplicitFunction / UnitSphere (autodiff, module chain)" "sanity-test-implicit.sh"
run_test "SphereMeshFixture (compute shader, binding layout)"      "sanity-test-sphere.sh"

echo ""
echo "════════════════════════════════════════"
echo " Summary: $PASS suite(s) passed, $FAIL failed"
echo "════════════════════════════════════════"
[[ $FAIL -eq 0 ]] && echo "All tests passed." || exit 1
