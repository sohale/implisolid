#!/bin/bash
# Sanity test: compile UnitSphere implicit function to SPIRV, WGSL, and HLSL.
# Verifies the module system, autodiff wiring, and access control are all intact.
#
# Usage:
#   bash scripts/sanity-test-implicit.sh
# Run from the repo root (implisolid/) or from slang/native/.
#
# What it tests:
#   - Module import chain: types → IImplicitFunction → UnitSphere
#   - __exported import (BoundingBox visible to callers)
#   - [Differentiable] + [NoDiffThis] + [BackwardDerivativeOf] on UnitSphere.eval
#   - gradient<T>() free function via bwd_diff(_evalForDiff<T>)
#   - Three compilation targets: spirv, wgsl, hlsl
#
# See docs/slang-learnings.md for explanation of every pattern used here.

set -euo pipefail

SLANGC=/dataneura/gpu-experimentations/experiments/20_slang_shaders/slang/build/RelWithDebInfo/bin/slangc

# Locate module directory relative to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="$SCRIPT_DIR/../implicit_function"
OUT_DIR="$(mktemp -d)"
trap 'rm -rf "$OUT_DIR"' EXIT

echo "slangc: $($SLANGC -v)"
echo "module dir: $MODULE_DIR"
echo "output dir: $OUT_DIR"
echo ""

# ---------------------------------------------------------------------------
# Inline test shader (no persistent test file in the source tree)
# ---------------------------------------------------------------------------
TEST_SHADER="$OUT_DIR/test_unit_sphere.slang"
cat > "$TEST_SHADER" << 'SLANG_EOF'
import UnitSphere;

struct TestOutput {
    float  f_inside;     // eval at p=(0.5,0,0) for r=1.5 sphere — expect > 0
    float3 grad_inside;  // gradient there — expect (-1, 0, 0)
    float  f_surface;    // eval at p=(1.5,0,0) — expect ~0
};

RWStructuredBuffer<TestOutput> gOutput;

[shader("compute")]
[numthreads(1, 1, 1)]
void main() {
    UnitSphere sphere = makeSphere(1.5f, float3(0.0f, 0.0f, 0.0f));

    float3 p_inside  = float3(0.5f, 0.0f, 0.0f);
    float3 p_surface = float3(1.5f, 0.0f, 0.0f);

    TestOutput out;
    out.f_inside    = sphere.eval(p_inside);
    out.grad_inside = gradient(sphere, p_inside);  // via autodiff in IImplicitFunction
    out.f_surface   = sphere.eval(p_surface);
    gOutput[0] = out;
}
SLANG_EOF

# ---------------------------------------------------------------------------
# Expected values (for documentation — not runtime-checked here):
#   f_inside  = 1.5² - 0.5² = 2.25 - 0.25 = 2.0         (> 0, inside ✓)
#   grad      = -2*(0.5-0, 0-0, 0-0) = (-1.0, 0.0, 0.0)  (inward ✓)
#   f_surface = 1.5² - 1.5² = 0.0                         (on surface ✓)
# ---------------------------------------------------------------------------

COMPILE="$SLANGC -I $MODULE_DIR $TEST_SHADER"

echo "--- SPIRV ---"
$COMPILE -target spirv -o "$OUT_DIR/out.spv"
echo "ok ($(wc -c < "$OUT_DIR/out.spv") bytes)"

echo ""
echo "--- WGSL ---"
$COMPILE -target wgsl -o "$OUT_DIR/out.wgsl"
echo "ok ($(wc -l < "$OUT_DIR/out.wgsl") lines)"

# Spot-check WGSL: verify hand-coded backward pass was used (not auto-generated)
if grep -q "UnitSphere_eval_bwd_0" "$OUT_DIR/out.wgsl"; then
    echo "autodiff wiring: hand-coded eval_bwd is present in WGSL output ✓"
else
    echo "WARNING: UnitSphere_eval_bwd_0 not found in WGSL — [BackwardDerivativeOf] may not have registered"
fi

echo ""
echo "--- HLSL ---"
$COMPILE -target hlsl -entry main -stage compute -o "$OUT_DIR/out.hlsl"
echo "ok ($(wc -l < "$OUT_DIR/out.hlsl") lines)"

echo ""
echo "All targets compiled successfully."
echo "See docs/slang-learnings.md for explanation of the patterns used."
