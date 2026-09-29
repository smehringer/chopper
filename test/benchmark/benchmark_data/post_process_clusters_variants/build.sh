#!/usr/bin/env bash
# Builds the `measure` binary. See README.md.
#
# Usage: CXX=<compiler> ./build.sh <output directory>
#
# Configures a chopper Release build in <output directory>/chopper, builds chopper_layout, applies variant_b.patch to
# a copy of src/layout/partition_user_bins.cpp and compiles measure.cpp with the flags that chopper uses for
# partition_user_bins.cpp.

set -Eeuo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: CXX=<compiler> $0 <output directory>" >&2
    exit 1
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CHOPPER_DIR=$(cd "${SCRIPT_DIR}/../../../.." && pwd)
OUT_DIR=$(mkdir -p "$1" && cd "$1" && pwd)
BUILD_DIR="${OUT_DIR}/chopper"

cmake -S "${CHOPPER_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build "${BUILD_DIR}" --target chopper_layout --parallel

PATCHED="${OUT_DIR}/partition_user_bins_variants.cpp"
patch --output="${PATCHED}" "${CHOPPER_DIR}/src/layout/partition_user_bins.cpp" < "${SCRIPT_DIR}/variant_b.patch"

# The compiler and flags of partition_user_bins.cpp, without -o <file>, -c and the source file.
mapfile -t COMMAND < <(python3 - "${BUILD_DIR}/compile_commands.json" <<'PYTHON'
import json, shlex, sys
entry = next(e for e in json.load(open(sys.argv[1])) if e["file"].endswith("src/layout/partition_user_bins.cpp"))
args = shlex.split(entry["command"])
result, skip = [], False
for arg in args:
    if skip:
        skip = False
    elif arg == "-o":
        skip = True
    elif arg != "-c" and arg != entry["file"]:
        result.append(arg)
print("\n".join(arg for arg in result if not arg.endswith("ccache")))
PYTHON
)

"${COMMAND[@]}" "${SCRIPT_DIR}/measure.cpp" "${PATCHED}" -o "${OUT_DIR}/measure" \
    "${BUILD_DIR}/lib/libchopper_layout.a" "${BUILD_DIR}/lib/libchopper_shared.a" "${BUILD_DIR}/lib/libhibf.a"

echo "Built ${OUT_DIR}/measure"
