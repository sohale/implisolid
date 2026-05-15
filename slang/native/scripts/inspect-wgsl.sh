#!/usr/bin/env bash
# Compile a Slang file to WGSL and display a structured summary of the output.
# Extracts: @binding declarations, struct definitions, workgroup sizes, entry points.
# Useful for verifying binding layout before writing a JS bridge.
#
# Usage:
#   inspect-wgsl.sh <file.slang> [entry]
#
# Examples:
#   inspect-wgsl.sh polygonisers/fixtures/sphere_surface_mesh_stub_for_test.slang genVertices
#   inspect-wgsl.sh polygonisers/fixtures/sphere_surface_mesh_stub_for_test.slang genFaces
#   inspect-wgsl.sh implicit_function/UnitSphere.slang

set -euo pipefail

SLANGC="${SLANGC:-/dataneura/gpu-experimentations/experiments/20_slang_shaders/slang/build/RelWithDebInfo/bin/slangc}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NATIVE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
IMPLICIT_DIR="$NATIVE_DIR/implicit_function"

if [[ $# -lt 1 ]]; then
    echo "Usage: $(basename "$0") <file.slang> [entry]" >&2
    exit 1
fi

FILE="$1"
ENTRY="${2:-}"

if [[ ! -f "$FILE" ]]; then
    if [[ -f "$NATIVE_DIR/$FILE" ]]; then
        FILE="$NATIVE_DIR/$FILE"
    else
        echo "error: file not found: $FILE" >&2
        exit 1
    fi
fi

WGSL="$(mktemp /tmp/inspect.XXXXXX.wgsl)"
trap "rm -f '$WGSL'" EXIT

CMD=("$SLANGC" -I "$IMPLICIT_DIR" "$FILE" -target wgsl -o "$WGSL")
[[ -n "$ENTRY" ]] && CMD+=(-entry "$ENTRY" -stage compute)

echo "Compiling: $FILE${ENTRY:+ (entry: $ENTRY)}"
"${CMD[@]}"
echo "Compiled: $(wc -l < "$WGSL") lines"
echo ""

echo "=== Bindings ==="
grep -n "@binding\|@group" "$WGSL" || echo "(none found)"

echo ""
echo "=== Structs ==="
grep -n "^struct " "$WGSL" || echo "(none found)"

echo ""
echo "=== Workgroup sizes ==="
grep -n "@workgroup_size\|workgroup_size" "$WGSL" || echo "(none found)"

echo ""
echo "=== Entry points (fn ...) ==="
grep -n "^fn " "$WGSL" || echo "(none found)"

echo ""
echo "=== Full WGSL output ==="
cat -n "$WGSL"
