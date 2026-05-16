#!/usr/bin/env bash
# Sanity tests for polygoniser fixtures and implicit-function modules.
# Compiles entry points (genVertices, genFaces, etc) to WGSL and verifies:
#   - Both entry points compile without errors
#   - Binding Layout: e.g., gParams=@binding(0), gVertices=@binding(1), gFaces=@binding(2), gNormals=@binding(3)
#   - Entry point names survive mangling: genVertices, genFaces
#   - Workgroup size: 8,8,1
#   - gVertices is array<f32>; gFaces is array<u32>
#
# Sections:
#   1. Sphere fixture    — genVertices / genFaces / genNormals
#   2. HalfSpace module  — eval, grad (constant), constructors, targets
#   3. CrispSubtract     — eval, piecewise grad, both inner grads in chain
#   4. Hemisphere fixture — genCapGeom / genDiskGeom / genCapFaces / genDiskFaces
#
# Usage:
#   bash scripts/sanity-test-sphere.sh

set -euo pipefail

SLANGC="${SLANGC:-/dataneura/gpu-experimentations/experiments/20_slang_shaders/slang/build/RelWithDebInfo/bin/slangc}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NATIVE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
IMPLICIT_DIR="$NATIVE_DIR/implicit_function"
FIXTURE_DIR="$NATIVE_DIR/polygonisers/fixtures"
SPHERE_FIXTURE="$FIXTURE_DIR/sphere_surface_mesh_stub_for_test.slang"
HEMI_FIXTURE="$FIXTURE_DIR/hemisphere_sharp_edge_fixture.slang"

TMPDIR_OWN="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_OWN"' EXIT

tmp() { mktemp "$TMPDIR_OWN/XXXXXX.$1"; }

PASS=0; FAIL=0

check() {
    local label="$1" result="$2"
    if [[ "$result" == "ok" ]]; then
        echo "  PASS  $label"; ((PASS++)) || true
    else
        echo "  FAIL  $label  ($result)"; ((FAIL++)) || true
    fi
}

compile_wgsl() {
    local src="$1" entry="$2" out="$3" extra="${4:-}"
    # shellcheck disable=SC2086
    "$SLANGC" -I "$IMPLICIT_DIR" $extra "$src" -target wgsl -entry "$entry" -stage compute -o "$out" 2>&1
}

compile_spirv() {
    local src="$1" out="$2" extra="${3:-}"
    # shellcheck disable=SC2086
    "$SLANGC" -I "$IMPLICIT_DIR" $extra "$src" -target spirv -o "$out" 2>&1
}

echo "includes: $IMPLICIT_DIR"
echo ""

# =============================================================================
# 1. SPHERE FIXTURE
# =============================================================================
echo "════════════════════════════════════════"
echo " 1. Sphere fixture"
echo "════════════════════════════════════════"
echo "fixture: $SPHERE_FIXTURE"

OUT_SV="$(tmp wgsl)"; OUT_SF="$(tmp wgsl)"; OUT_SN="$(tmp wgsl)"
SCOMPILE="$SLANGC -I $IMPLICIT_DIR $SPHERE_FIXTURE -target wgsl -stage compute"

echo ""
echo "Compiling..."
$SCOMPILE -entry genVertices -o "$OUT_SV" && check "genVertices compiles" "ok" || check "genVertices compiles" "slangc failed"
$SCOMPILE -entry genFaces    -o "$OUT_SF" && check "genFaces compiles"    "ok" || check "genFaces compiles"    "slangc failed"
$SCOMPILE -entry genNormals  -o "$OUT_SN" && check "genNormals compiles"  "ok" || check "genNormals compiles"  "slangc failed"

echo ""
echo "--- genVertices ---"
check "binding(0) params"           "$(grep -q '@binding(0)' "$OUT_SV"                       && echo ok || echo 'not found')"
check "binding(1) vertices array<f32>" "$(grep -q '@binding(1)' "$OUT_SV" && grep -q 'array<f32>' "$OUT_SV" && echo ok || echo 'not found')"
check "fn genVertices"              "$(grep -q 'fn genVertices'    "$OUT_SV"                  && echo ok || echo 'not found')"
check "workgroup_size(8, 8, 1)"     "$(grep -q 'workgroup_size(8, 8, 1)' "$OUT_SV"           && echo ok || echo 'not found')"
check "no binding(2) or (3)"        "$(! grep -q '@binding(2)\|@binding(3)' "$OUT_SV"        && echo ok || echo 'unexpectedly present')"

echo ""
echo "--- genFaces ---"
check "binding(0) params"           "$(grep -q '@binding(0)'  "$OUT_SF"                      && echo ok || echo 'not found')"
check "binding(2) faces array<u32>" "$(grep -q '@binding(2)'  "$OUT_SF" && grep -q 'array<u32>' "$OUT_SF" && echo ok || echo 'not found')"
check "fn genFaces"                 "$(grep -q 'fn genFaces'   "$OUT_SF"                     && echo ok || echo 'not found')"
check "no binding(1) or (3)"        "$(! grep -q '@binding(1)\|@binding(3)' "$OUT_SF"        && echo ok || echo 'unexpectedly present')"

echo ""
echo "--- genNormals ---"
check "binding(0) params"           "$(grep -q '@binding(0)'  "$OUT_SN"                      && echo ok || echo 'not found')"
check "binding(1) vertices"         "$(grep -q '@binding(1)'  "$OUT_SN"                      && echo ok || echo 'not found')"
check "binding(3) normals"          "$(grep -q '@binding(3)'  "$OUT_SN"                      && echo ok || echo 'not found')"
check "fn genNormals"               "$(grep -q 'fn genNormals' "$OUT_SN"                     && echo ok || echo 'not found')"
check "workgroup_size(8, 8, 1)"     "$(grep -q 'workgroup_size(8, 8, 1)' "$OUT_SN"           && echo ok || echo 'not found')"
check "no binding(2)"               "$(! grep -q '@binding(2)' "$OUT_SN"                     && echo ok || echo 'unexpectedly present')"
check "autodiff: UnitSphere_eval_bwd" "$(grep -q 'UnitSphere_eval_bwd_0' "$OUT_SN"           && echo ok || echo 'not found')"
check "autodiff: gradient fn"       "$(grep -q 'fn gradient_0' "$OUT_SN"                     && echo ok || echo 'not found')"
check "outward normal (negate grad)" "$(grep -q 'vec3<f32>(0) - gradient_0\|-.*gradient' "$OUT_SN" && echo ok || echo 'not found')"

# =============================================================================
# 2. HALFSPACE MODULE
# =============================================================================
echo ""
echo "════════════════════════════════════════"
echo " 2. HalfSpace module"
echo "════════════════════════════════════════"

# Inline test shader
HS_TEST="$(tmp slang)"
cat > "$HS_TEST" << 'SLANG_EOF'
import HalfSpace;
RWStructuredBuffer<float> gOut;
[shader("compute")][numthreads(1,1,1)]
void main() {
    // z-positive half-space: f(p) = p.z, gradient = (0,0,1)
    HalfSpace zpos = makeZPositiveHalfSpace();
    float3 p = float3(1.0f, 2.0f, 3.0f);
    gOut[0] = zpos.eval(p);          // = 3.0  (p.z)
    float3 g = gradient(zpos, p);    // = (0,0,1)  — constant, linear function
    gOut[1] = g.x; gOut[2] = g.y; gOut[3] = g.z;

    // z-negative half-space: f(p) = -p.z, gradient = (0,0,-1)
    HalfSpace zneg = makeZNegativeHalfSpace();
    gOut[4] = zneg.eval(p);          // = -3.0
    float3 gn = gradient(zneg, p);   // = (0,0,-1)
    gOut[5] = gn.x; gOut[6] = gn.y; gOut[7] = gn.z;
}
SLANG_EOF

OUT_HS_WGSL="$(tmp wgsl)"
OUT_HS_SPV="$(tmp spv)"

echo ""
echo "Compiling HalfSpace test..."
compile_wgsl "$HS_TEST" main "$OUT_HS_WGSL" && check "HalfSpace compiles to WGSL"  "ok" || check "HalfSpace compiles to WGSL"  "slangc failed"
compile_spirv "$HS_TEST" "$OUT_HS_SPV"       && check "HalfSpace compiles to SPIRV" "ok" || check "HalfSpace compiles to SPIRV" "slangc failed"

echo ""
echo "--- HalfSpace WGSL checks ---"
check "HalfSpace_eval_bwd_0 (hand-coded grad)" "$(grep -q 'HalfSpace_eval_bwd_0'         "$OUT_HS_WGSL" && echo ok || echo 'not found')"
check "makeZPositiveHalfSpace_0 constructor"    "$(grep -q 'makeZPositiveHalfSpace_0'     "$OUT_HS_WGSL" && echo ok || echo 'not found')"
check "makeZNegativeHalfSpace_0 constructor"    "$(grep -q 'makeZNegativeHalfSpace_0'     "$OUT_HS_WGSL" && echo ok || echo 'not found')"
check "gradient fn present"                     "$(grep -q 'fn gradient_0'                "$OUT_HS_WGSL" && echo ok || echo 'not found')"
# Gradient of linear function is constant (normal field) — no conditional branch in eval_bwd.
# WGSL mangles 'normal' → 'normal_0'; function boundary found via next '\nfn '.
check "grad is constant (no if in eval_bwd)"    "$(python3 -c "
import re
txt = open('$OUT_HS_WGSL').read()
idx = txt.find('fn HalfSpace_eval_bwd_0')
if idx == -1:
    print('HalfSpace_eval_bwd_0 not found')
else:
    rest  = txt[idx:]
    end   = re.search(r'\nfn ', rest[10:])
    body  = rest[: end.start() + 10] if end else rest
    has_normal = 'normal_0' in body
    has_if     = bool(re.search(r'\bif\s*\(', body))
    print('ok' if has_normal and not has_if else f'conditional={has_if} normal_field={has_normal}')
" 2>/dev/null || echo 'parse-error')"

# =============================================================================
# 3. CRISPSUBTRACT MODULE
# =============================================================================
echo ""
echo "════════════════════════════════════════"
echo " 3. CrispSubtract<UnitSphere, HalfSpace>"
echo "════════════════════════════════════════"

CS_TEST="$(tmp slang)"
cat > "$CS_TEST" << 'SLANG_EOF'
import UnitSphere;
import HalfSpace;
import CrispSubtract;
RWStructuredBuffer<float> gOut;
[shader("compute")][numthreads(1,1,1)]
void main() {
    // hemisphere = sphere - negative-Z space
    UnitSphere sphere = makeSphere(1.0f, float3(0.0f));
    HalfSpace  negz   = makeZNegativeHalfSpace();
    var hemi = makeCrispSubtract(sphere, negz);

    float3 p_cap  = float3(0.0f, 0.0f, 0.9f);  // on curved cap, near north pole
    float3 p_disk = float3(0.9f, 0.0f, 0.0f);  // on equator, disk side
    gOut[0] = hemi.eval(p_cap);
    gOut[1] = hemi.eval(p_disk);
    float3 g_cap  = gradient(hemi, p_cap);      // should be sphere-gradient-like
    float3 g_disk = gradient(hemi, p_disk);     // should be (0,0,-1) direction
    gOut[2] = g_cap.x;  gOut[3] = g_cap.y;  gOut[4] = g_cap.z;
    gOut[5] = g_disk.x; gOut[6] = g_disk.y; gOut[7] = g_disk.z;
}
SLANG_EOF

OUT_CS_WGSL="$(tmp wgsl)"
OUT_CS_SPV="$(tmp spv)"

echo ""
echo "Compiling CrispSubtract test..."
compile_wgsl "$CS_TEST" main "$OUT_CS_WGSL" && check "CrispSubtract compiles to WGSL"  "ok" || check "CrispSubtract compiles to WGSL"  "slangc failed"
compile_spirv "$CS_TEST" "$OUT_CS_SPV"       && check "CrispSubtract compiles to SPIRV" "ok" || check "CrispSubtract compiles to SPIRV" "slangc failed"

echo ""
echo "--- CrispSubtract WGSL checks ---"
check "CrispSubtract_eval_bwd_0 (manual grad)" "$(grep -q 'CrispSubtract_eval_bwd_0'     "$OUT_CS_WGSL" && echo ok || echo 'not found')"
check "UnitSphere_eval_bwd_0 in chain"          "$(grep -q 'UnitSphere_eval_bwd_0'         "$OUT_CS_WGSL" && echo ok || echo 'not found')"
check "HalfSpace_eval_bwd_0 in chain"           "$(grep -q 'HalfSpace_eval_bwd_0'          "$OUT_CS_WGSL" && echo ok || echo 'not found')"
check "min() in eval output"                    "$(grep -q '\bmin(' "$OUT_CS_WGSL"                        && echo ok || echo 'not found')"
check "makeCrispSubtract_0 constructor"         "$(grep -q 'makeCrispSubtract_0'           "$OUT_CS_WGSL" && echo ok || echo 'not found')"
# Piecewise gradient — CrispSubtract.eval_bwd must contain an 'if' (the branch on f_a < -f_b).
# Function boundary: next '\nfn ' after the start. WGSL emits 'if((' not 'if('.
check "piecewise grad: if branch in eval_bwd"   "$(python3 -c "
import re
txt = open('$OUT_CS_WGSL').read()
idx = txt.find('fn CrispSubtract_eval_bwd_0')
if idx == -1:
    print('CrispSubtract_eval_bwd_0 not found')
else:
    rest = txt[idx:]
    end  = re.search(r'\nfn ', rest[10:])
    body = rest[: end.start() + 10] if end else rest
    print('ok' if bool(re.search(r'\bif\s*\(', body)) else 'no conditional found')
" 2>/dev/null || echo 'parse-error')"
# Both inner gradients (gradient_0 = UnitSphere, gradient_1 = HalfSpace) must appear in eval_bwd.
check "both inner gradients called from eval_bwd" "$(python3 -c "
import re
txt = open('$OUT_CS_WGSL').read()
idx = txt.find('fn CrispSubtract_eval_bwd_0')
if idx == -1:
    print('CrispSubtract_eval_bwd_0 not found')
else:
    rest = txt[idx:]
    end  = re.search(r'\nfn ', rest[10:])
    body = rest[: end.start() + 10] if end else rest
    g0 = 'gradient_0' in body
    g1 = 'gradient_1' in body
    print('ok' if g0 and g1 else f'missing: gradient_0={g0} gradient_1={g1}')
" 2>/dev/null || echo 'parse-error')"
# The gradient of -f_B branch must negate the HalfSpace gradient (verified via minus sign).
check "negation of inner grad (-gradient_b)" "$(grep -q '\-.*gradient\|gradient.*neg\|vec3<f32>(0\.0' "$OUT_CS_WGSL" && echo ok || echo 'not found')"

# =============================================================================
# 4. HEMISPHERE FIXTURE
# =============================================================================
echo ""
echo "════════════════════════════════════════"
echo " 4. Hemisphere fixture"
echo "════════════════════════════════════════"
echo "fixture: $HEMI_FIXTURE"

OUT_HCG="$(tmp wgsl)"; OUT_HDG="$(tmp wgsl)"
OUT_HCF="$(tmp wgsl)"; OUT_HDF="$(tmp wgsl)"
HCOMPILE="$SLANGC -I $IMPLICIT_DIR $HEMI_FIXTURE -target wgsl -stage compute"

echo ""
echo "Compiling..."
$HCOMPILE -entry genCapGeom   -o "$OUT_HCG" && check "genCapGeom compiles"   "ok" || check "genCapGeom compiles"   "slangc failed"
$HCOMPILE -entry genDiskGeom  -o "$OUT_HDG" && check "genDiskGeom compiles"  "ok" || check "genDiskGeom compiles"  "slangc failed"
$HCOMPILE -entry genCapFaces  -o "$OUT_HCF" && check "genCapFaces compiles"  "ok" || check "genCapFaces compiles"  "slangc failed"
$HCOMPILE -entry genDiskFaces -o "$OUT_HDF" && check "genDiskFaces compiles" "ok" || check "genDiskFaces compiles" "slangc failed"

echo ""
echo "--- Binding layout (geom={0,1,3}, faces={0,2}) ---"
for f in "$OUT_HCG" "$OUT_HDG"; do
    name=$(grep -o 'fn gen[A-Za-z]*' "$f" | head -1)
    check "$name: binding(0) params"        "$(grep -q '@binding(0)' "$f"                        && echo ok || echo 'not found')"
    check "$name: binding(1) vertices"      "$(grep -q '@binding(1)' "$f"                        && echo ok || echo 'not found')"
    check "$name: binding(3) normals"       "$(grep -q '@binding(3)' "$f"                        && echo ok || echo 'not found')"
    check "$name: no binding(2) (faces)"    "$(! grep -q '@binding(2)' "$f"                      && echo ok || echo 'unexpectedly present')"
done
for f in "$OUT_HCF" "$OUT_HDF"; do
    name=$(grep -o 'fn gen[A-Za-z]*' "$f" | head -1)
    check "$name: binding(0) params"        "$(grep -q '@binding(0)' "$f"                        && echo ok || echo 'not found')"
    check "$name: binding(2) faces"         "$(grep -q '@binding(2)' "$f"                        && echo ok || echo 'not found')"
    check "$name: no binding(1) or (3)"     "$(! grep -q '@binding(1)\|@binding(3)' "$f"         && echo ok || echo 'unexpectedly present')"
done

echo ""
echo "--- Workgroup sizes ---"
check "genCapGeom:   workgroup(8,8,1)"   "$(grep -q 'workgroup_size(8, 8, 1)'   "$OUT_HCG" && echo ok || echo 'not found')"
check "genDiskGeom:  workgroup(64,1,1)"  "$(grep -q 'workgroup_size(64, 1, 1)'  "$OUT_HDG" && echo ok || echo 'not found')"
check "genCapFaces:  workgroup(8,8,1)"   "$(grep -q 'workgroup_size(8, 8, 1)'   "$OUT_HCF" && echo ok || echo 'not found')"
check "genDiskFaces: workgroup(64,1,1)"  "$(grep -q 'workgroup_size(64, 1, 1)'  "$OUT_HDF" && echo ok || echo 'not found')"

echo ""
echo "--- Gradient correctness per section ---"
check "genCapGeom:  UnitSphere gradient (sphere normals)" \
    "$(grep -q 'UnitSphere_eval_bwd_0' "$OUT_HCG"                                         && echo ok || echo 'not found')"
check "genCapGeom:  outward sphere normal (negate grad)" \
    "$(grep -q 'vec3<f32>(0) - gradient\|-.*gradient' "$OUT_HCG"                          && echo ok || echo 'not found')"
check "genDiskGeom: HalfSpace gradient (disk normals)" \
    "$(grep -q 'HalfSpace_eval_bwd_0' "$OUT_HDG"                                          && echo ok || echo 'not found')"
check "genDiskGeom: outward disk normal (negate z-up grad → 0,0,-1)" \
    "$(grep -q 'vec3<f32>(0) - gradient\|-.*gradient' "$OUT_HDG"                          && echo ok || echo 'not found')"
check "genCapGeom:  NOT HalfSpace (cap uses sphere only)" \
    "$(! grep -q 'HalfSpace_eval_bwd_0' "$OUT_HCG"                                        && echo ok || echo 'HalfSpace unexpectedly present')"
check "genDiskGeom: NOT UnitSphere (disk uses halfspace only)" \
    "$(! grep -q 'UnitSphere_eval_bwd_0' "$OUT_HDG"                                       && echo ok || echo 'UnitSphere unexpectedly present')"

# =============================================================================
echo ""
echo "════════════════════════════════════════"
echo " Results: $PASS passed, $FAIL failed"
echo "════════════════════════════════════════"
[[ $FAIL -eq 0 ]] && echo "All checks passed." || { echo "Some checks failed."; exit 1; }
