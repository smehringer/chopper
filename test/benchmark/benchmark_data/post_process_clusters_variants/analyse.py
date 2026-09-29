#!/usr/bin/env python3
# ---------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
# ---------------------------------------------------------------------------------------------------
"""Summarises the output of run.sh. See README.md.

Usage: analyse.py results_small.tsv [results_large.tsv ...]

For each file, prints the ratio variant B / variant A per input and, over the inputs whose top-level partitions differ,
the number of inputs where B is better/worse, the geometric mean and range of B/A, a two-sided sign test and, if SciPy
is installed, a two-sided Wilcoxon signed-rank test on log(B/A). Lower is better for all metrics.
"""

import csv
import math
import sys

try:
    from scipy.stats import wilcoxon
except ImportError:
    wilcoxon = None

METRICS = [("topmax", "largest top-level technical bin"), ("size", "expected HIBF size"), ("cost", "expected query cost")]


def sign_test(better: int, worse: int) -> float:
    """Exact two-sided sign test."""
    n = better + worse
    if n == 0:
        return 1.0
    tail = sum(math.comb(n, i) for i in range(min(better, worse) + 1)) / 2**n
    return min(1.0, 2 * tail)


def summarise(rows: list[dict], title: str) -> None:
    print(f"## {title}\n")
    print("| n | families | sigma | seed | tmax | s | same | topmax B/A | size B/A | cost B/A |")
    print("|--:|--:|--:|--:|--:|--:|:-:|--:|--:|--:|")
    for r in rows:
        ratios = " | ".join(f"{float(r[m + '_B']) / float(r[m + '_A']):.4f}" for m, _ in METRICS)
        same = "yes" if r["same_top"] == "1" else "no"
        print(f"| {r['n']} | {r['families']} | {r['sigma']} | {r['seed']} | {r['tmax']} | {r['s_A']} | {same} | {ratios} |")

    differing = [r for r in rows if r["same_top"] == "0"]
    print(f"\nInputs: {len(rows)}. Top-level partitions differ in {len(differing)}.\n")
    if not differing:
        return
    print("| Metric | B better | B worse | equal | geo-mean B/A | range B/A | sign test p | Wilcoxon p |")
    print("|---|--:|--:|--:|--:|--:|--:|--:|")
    for metric, name in METRICS:
        logs = [math.log(float(r[metric + "_B"]) / float(r[metric + "_A"])) for r in differing]
        nonzero = [x for x in logs if abs(x) > 1e-12]
        better = sum(x < 0 for x in nonzero)
        worse = sum(x > 0 for x in nonzero)
        geo_mean = math.exp(sum(logs) / len(logs))
        wilcoxon_p = f"{wilcoxon(nonzero).pvalue:.2f}" if wilcoxon and nonzero else "n/a"
        print(
            f"| {name} | {better} | {worse} | {len(logs) - len(nonzero)} | {geo_mean:.4f} | "
            f"{math.exp(min(logs)):.3f} - {math.exp(max(logs)):.3f} | {sign_test(better, worse):.2f} | {wilcoxon_p} |"
        )
    print()


def main() -> None:
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for path in sys.argv[1:]:
        with open(path, newline="") as f:
            summarise(list(csv.DictReader(f, delimiter="\t")), path)


if __name__ == "__main__":
    main()
