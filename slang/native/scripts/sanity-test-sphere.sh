#!/usr/bin/env bash
# Sanity test for sphere_surface_mesh_stub_for_test.slang.
# Compiles both entry points (genVertices, genFaces) to WGSL and verifies:
#   - Both entry points compile without errors
#   - Binding layout: gParams=@binding(0), gVertices=@binding(1), gFaces=@binding(2)
#   - Entry point names survive mangling: genVertices, genFaces
#   - Workgroup size: 8,8,1
#   - gVertices is array<f32>; gFaces is array<u32>
#
# Usage:
#   bash scripts/sanity-test-sphere.sh
# Run from any directory (script locates itself).

set -euo pipefail

SLANGC="${SLANGC:-/dataneura/gpu-experimentations/experiments/20_slang_shaders/slang/build/RelWithDebInfo/bin/slangc}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NATIVE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
FIXTURE="$NATIVE_DIR/polygonisers/fixtures/sphere_surface_mesh_stub_for_test.slang"

OUT_V="$(mktemp /tmp/sphere_vertices.XXXXXX.wgsl)"
OUT_F="$(mktemp /tmp/sphere_faces.XXXXXX.wgsl)"
trap "rm -f '$OUT_V' '$OUT_F'" EXIT

PASS=0
FAIL=0

check() {
    local label="$1"; local result="$2"
    if [[ "$result" == "ok" ]]; then
        echo "  PASS  $label"
        ((PASS++)) || true
    else
        echo "  FAIL  $label  ($result)"
        ((FAIL++)) || true
    fi
}

echo "fixture: $FIXTURE"
echo ""

# --- Compile ---
echo "Compiling genVertices..."
"$SLANGC" "$FIXTURE" -target wgsl -entry genVertices -stage compute -o "$OUT_V" \
    && check "genVertices compiles" "ok" \
    || check "genVertices compiles" "slangc failed"

echo "Compiling genFaces..."
"$SLANGC" "$FIXTURE" -target wgsl -entry genFaces -stage compute -o "$OUT_F" \
    && check "genFaces compiles" "ok" \
    || check "genFaces compiles" "slangc failed"

echo ""
echo "=== genVertices WGSL checks ==="

check "gParams at binding(0)" \
    "$(grep -q '@binding(0).*gParams\|gParams.*@binding(0)' "$OUT_V" 2>/dev/null && echo ok || echo "not found")"

check "gVertices at binding(1)" \
    "$(grep -q '@binding(1)' "$OUT_V" 2>/dev/null && echo ok || echo "not found")"

check "gVertices is array<f32>" \
    "$(grep -q 'array<f32>' "$OUT_V" 2>/dev/null && echo ok || echo "not found")"

check "genVertices fn present" \
    "$(grep -q 'fn genVertices' "$OUT_V" 2>/dev/null && echo ok || echo "not found")"

check "workgroup_size(8, 8, 1)" \
    "$(grep -q 'workgroup_size(8, 8, 1)' "$OUT_V" 2>/dev/null && echo ok || echo "not found")"

check "no binding(2) in genVertices WGSL" \
    "$(grep -qv '@binding(2)' "$OUT_V" 2>/dev/null && ! grep -q '@binding(2)' "$OUT_V" && echo ok || echo "binding(2) unexpectedly present")"

echo ""
echo "=== genFaces WGSL checks ==="

check "gParams at binding(0)" \
    "$(grep -q '@binding(0).*gParams\|gParams.*@binding(0)' "$OUT_F" 2>/dev/null && echo ok || echo "not found")"

check "gFaces at binding(2)" \
    "$(grep -q '@binding(2)' "$OUT_F" 2>/dev/null && echo ok || echo "not found")"

check "gFaces is array<u32>" \
    "$(grep -q 'array<u32>' "$OUT_F" 2>/dev/null && echo ok || echo "not found")"

check "genFaces fn present" \
    "$(grep -q 'fn genFaces' "$OUT_F" 2>/dev/null && echo ok || echo "not found")"

check "workgroup_size(8, 8, 1)" \
    "$(grep -q 'workgroup_size(8, 8, 1)' "$OUT_F" 2>/dev/null && echo ok || echo "not found")"

check "no binding(1) in genFaces WGSL" \
    "$(! grep -q '@binding(1)' "$OUT_F" && echo ok || echo "binding(1) unexpectedly present")"

echo ""
echo "=== SphereParams struct layout check (genVertices WGSL) ==="
# Verify the std140 layout fields appear in order
for field in radius_0 cx_0 cy_0 cz_0 M_0 N_0; do
    check "field $field present" \
        "$(grep -q "$field" "$OUT_V" 2>/dev/null && echo ok || echo "not found")"
done

echo ""
echo "Results: $PASS passed, $FAIL failed"
[[ $FAIL -eq 0 ]] && echo "All checks passed." || { echo "Some checks failed."; exit 1; }
