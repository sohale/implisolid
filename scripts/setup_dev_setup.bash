#!/bin/bash
set -euo pipefail
set -x

# Dev Setup: Sets up & checks the dev evironment for: vscode, clangd LSP, etc.
# run check_dev_setup.bash to generate compile_commands.json in build/
# Disable Microsoft IntelliSense diagnostics: Id: ms-vscode.cpptools , Id: ms-vscode.cpptools-extension-pack.
# and let clangd be the only authority.



export ORIG_REPO_ROOT="$(git rev-parse --show-toplevel)"
export E2E="$ORIG_REPO_ROOT/e2e-sandbox-temp"
echo "Repo root        : $ORIG_REPO_ROOT"


CLANGD="/usr/bin/clangd-18"
CXX="clang++-18"
# Build-time only (not relied upon by clangd)
export CC="clang-18"
export CXX


BOOST="$ORIG_REPO_ROOT/build/lib/boost_1_75_0"
EIGEN="$ORIG_REPO_ROOT/build/lib/eigen"




SRC1="$ORIG_REPO_ROOT/js_iteration_1/mcc2.cpp"
test -f "$SRC1"

COMPILE_DB="$ORIG_REPO_ROOT/build/compile_commands.json"

cat > "$COMPILE_DB" <<EOF
[
  {
    "directory": "$ORIG_REPO_ROOT",
    "file": "$SRC1",
    "command": "$CXX -std=c++17 -I$BOOST -I$EIGEN js_iteration_1/mcc2.cpp"
  }
]
EOF

echo "compile_commands.json written to:"
echo "  $COMPILE_DB"
echo "Contents: ======="
batcat -pp "$COMPILE_DB" || cat -pp "$COMPILE_DB"
echo -e "=======\n"


"$CLANGD" --version
"$CXX" --version
"$CLANGD" --version

"$CLANGD" --check=$SRC1 --compile-commands-dir=build

cd "$ORIG_REPO_ROOT"
# test -f .clangd # not anymore
test -f ".vscode/settings.json"
test -f "$COMPILE_DB"
