#!/usr/bin/env bash
# General-purpose Slang compiler wrapper.
#
# Usage:
#   compile.sh <file.slang> [target] [entry] [output]
#
#   target  — wgsl (default), spirv, hlsl, cuda, cpp
#   entry   — compute entry point name (required for hlsl; optional for wgsl/spirv)
#   output  — path for compiled output (default: /tmp/out.<target>)
#
# Examples:
#   compile.sh implicit_function/UnitSphere.slang
#   compile.sh implicit_function/UnitSphere.slang wgsl
#   compile.sh implicit_function/UnitSphere.slang hlsl main /tmp/sphere.hlsl
#   compile.sh polygonisers/fixtures/sphere_surface_mesh_stub_for_test.slang wgsl genVertices
#
# The -I flag is always set to slang/native/implicit_function/ for module resolution.

set -euo pipefail

SLANGC="${SLANGC:-/dataneura/gpu-experimentations/experiments/20_slang_shaders/slang/build/RelWithDebInfo/bin/slangc}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NATIVE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
IMPLICIT_DIR="$NATIVE_DIR/implicit_function"

if [[ $# -lt 1 ]]; then
    echo "Usage: $(basename "$0") <file.slang> [target] [entry] [output]" >&2
    exit 1
fi

FILE="$1"
TARGET="${2:-wgsl}"
ENTRY="${3:-}"
OUTPUT="${4:-/tmp/out.$TARGET}"

# Resolve relative paths from either repo root or native dir
if [[ ! -f "$FILE" ]]; then
    if [[ -f "$NATIVE_DIR/$FILE" ]]; then
        FILE="$NATIVE_DIR/$FILE"
    else
        echo "error: file not found: $FILE" >&2
        exit 1
    fi
fi

CMD=("$SLANGC" -I "$IMPLICIT_DIR" "$FILE" -target "$TARGET" -o "$OUTPUT")
[[ -n "$ENTRY" ]] && CMD+=(-entry "$ENTRY" -stage compute)

echo "slangc $(${SLANGC} -v 2>&1 | head -1)"
echo "target:  $TARGET"
echo "entry:   ${ENTRY:-(auto)}"
echo "output:  $OUTPUT"
echo "cmd:     ${CMD[*]}"
echo ""

"${CMD[@]}"
echo "ok — $(wc -c < "$OUTPUT") bytes written to $OUTPUT"
