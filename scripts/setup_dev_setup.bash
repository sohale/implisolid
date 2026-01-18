#!/bin/bash
set -euo pipefail
set -x

# Dev Setup: Sets up & checks the dev evironment for: vscode, clangd LSP, etc.
# run check_dev_setup.bash to generate compile_commands.json in build/
# Disable Microsoft IntelliSense diagnostics: Id: ms-vscode.cpptools , Id: ms-vscode.cpptools-extension-pack.
# and let clangd be the only authority.

# llvm-vs-code-extensions.vscode-clangd

# Note that this does not use the e2e folder that keeps a copy (scripted build), but the actual original repo's files (as required for `clangd` LSP to work in a useful way).
export ORIG_REPO_ROOT="$(git rev-parse --show-toplevel)"
# export E2E="$ORIG_REPO_ROOT/e2e-sandbox-temp"
echo "Repo root        : $ORIG_REPO_ROOT" # will be: "/dataneura/3d/implisolid"

# LLVERSION="18"
LLVERSION="21"

CLANGD="/usr/bin/clangd-${LLVERSION}"
CXX="clang++-${LLVERSION}"
# Build-time only (not relied upon by clangd)
export CC="clang-${LLVERSION}"
export CXX


BOOST="$ORIG_REPO_ROOT/build/lib/boost_1_75_0"
EIGEN="$ORIG_REPO_ROOT/build/lib/eigen"

# Emscripten:
EMSDK_SURROGATE="build/emscripten_emsdk_surrogate_for_clangd"
export EMSDK_VERSION="3.1.14"
# todo: unify with another similat defnition

mkdir -p "$EMSDK_SURROGATE"
docker run --rm \
  -u 1000:1000 \
  -v "${ORIG_REPO_ROOT}/${EMSDK_SURROGATE}:/out" \
  emscripten/emsdk:${EMSDK_VERSION} \
  bash -c '
    set -euo pipefail
    cp -a /emsdk/upstream/emscripten/cache/sysroot/include /out/
    echo "I am in emscripten"
    # find /emsdk
  '
# verify the purpose is fulfilled:
test -f "${EMSDK_SURROGATE}/include/emscripten.h"
test -f "${EMSDK_SURROGATE}/include/c++/v1/vector"

# Wrong: /emsdk/upstream/emscripten/system/include
# Right: /emsdk/upstream/emscripten/cache/sysroot/include
# Reason: Recent versions of the SDK (3.x) moved to a model where only the headers under the cache/sysroot include tree are intended for compilation. It contains sanitized copies of standard headers (libc/libc++), platform headers, and Emscripten API headers that actually work with both native and cross-compilation modes.

# A dummy file will be required: clangd_basic_data_structures.cpp

CLANGD_TU_DIR="$ORIG_REPO_ROOT/build/clangd_TUs"
CLANG_TU1="$CLANGD_TU_DIR/clangd_basic_data_structures.cpp"
mkdir -p "$CLANGD_TU_DIR"
echo '
// Dummy file for clangd errors.
#include <vector>
#include <cassert>
#include "js_iteration_2/basic_data_structures.hpp"
#include "js_iteration_2/object_factory.hpp"
' > "${CLANG_TU1}"

SRC1="$ORIG_REPO_ROOT/js_iteration_1/mcc2.cpp"
test -f "$SRC1"

COMPILE_DB="$ORIG_REPO_ROOT/build/compile_commands.json"

CLANGD_ARGS=(
  -std=c++17
  -I$BOOST
  -I$EIGEN
  --sysroot=${ORIG_REPO_ROOT}/${EMSDK_SURROGATE}
  -isystem ${ORIG_REPO_ROOT}/${EMSDK_SURROGATE}/include/c++/v1
  --target=wasm32-unknown-emscripten
  -nostdinc++
)

cat > "$COMPILE_DB" <<EOF
[
  {
    "directory": "$ORIG_REPO_ROOT",
    "file": "$SRC1",
    "command": "$CXX  ${CLANGD_ARGS[*]}  js_iteration_1/mcc2.cpp"
  },
  {
    "directory": "$ORIG_REPO_ROOT",
    "file": "$CLANG_TU1",
    "command": "$CXX  ${CLANGD_ARGS[*]}   ${CLANG_TU1}"
  }
]
EOF

# was:
#. "$CXX -std=c++17  -I$BOOST -I$EIGEN  -isystem ${ORIG_REPO_ROOT}/${EMSDK_SURROGATE}/include -isystem ${ORIG_REPO_ROOT}/${EMSDK_SURROGATE}/include/c++/v1   --target=wasm32-unknown-emscripten -nostdinc++ js_iteration_1/mcc2.cpp"
# This works cleanly: 0 errors:
# "$CXX -std=c++17  -I$BOOST -I$EIGEN  --sysroot=${ORIG_REPO_ROOT}/${EMSDK_SURROGATE} -isystem ${ORIG_REPO_ROOT}/${EMSDK_SURROGATE}/include/c++/v1   --target=wasm32-unknown-emscripten -nostdinc++ js_iteration_1/mcc2.cpp"

echo "compile_commands.json written to:"
echo "  $COMPILE_DB"
echo "Contents: ======="
batcat -pp "$COMPILE_DB" || cat -pp "$COMPILE_DB"
echo -e "=======\n"


"$CLANGD" --version
# "$CXX" --version

"$CLANGD" --check=$SRC1 --compile-commands-dir=build

cd "$ORIG_REPO_ROOT"
# test -f .clangd # not anymore
test -f ".vscode/settings.json"
test -f "$COMPILE_DB"


# Now you can use, even in commandline, :
# /usr/bin/clangd-21  --compile-commands-dir=build --check=js_iteration_2/implicit_function/cube.hpp
