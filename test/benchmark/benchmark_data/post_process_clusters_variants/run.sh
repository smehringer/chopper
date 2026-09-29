#!/usr/bin/env bash
# Runs both input grids of README.md and writes results_small.tsv and results_large.tsv to the current directory.
#
# Usage: ./run.sh <path to measure>

set -Eeuo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <path to measure>" >&2
    exit 1
fi

MEASURE=$1
HEADER="n\tfamilies\tsigma\tseed\tmedian\tmin\tmax\ttmax\ts_A\ts_B\tsame_top\ttopmax_A\ttopmax_B\tsize_A\tsize_B\tcost_A\tcost_B\tdepth_A\tdepth_B\tt_sketch\tt_A\tt_B"

# Small grid: 400 and 1000 user bins, median 8000 k-mers, 3 seeds each.
printf '%b\n' "${HEADER}" > results_small.tsv
for n in 400 1000; do
    for family_divisor in 1 4 25; do
        for sigma in 0.7 1.4; do
            for seed in 1 2 3; do
                "${MEASURE}" "${n}" "${family_divisor}" "${sigma}" "${seed}" 8000 5000 400000 8 >> results_small.tsv
            done
        done
    done
done

# Large grid: 10000, 25000 and 50000 small user bins, median 4000 k-mers, 3, 2 and 1 seeds.
printf '%b\n' "${HEADER}" > results_large.tsv
for spec in "10000 3" "25000 2" "50000 1"; do
    read -r n seeds <<< "${spec}"
    for family_divisor in 5 50; do
        for sigma in 1.4 2.0; do
            for seed in $(seq 1 "${seeds}"); do
                "${MEASURE}" "${n}" "${family_divisor}" "${sigma}" "${seed}" 4000 2500 2000000 32 >> results_large.tsv
            done
        done
    done
done
