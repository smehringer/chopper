## Seeding clusters in `post_process_clusters`

Measurement behind the decision to keep how `post_process_clusters` (`src/layout/partition_user_bins.cpp`) orders the
clusters for the fast layout.

### Question

`lsh_sim_approach` seeds the merged technical bins with the first clusters, then assigns all other clusters greedily by
similarity (`find_best_partition`). Before that, `post_process_clusters` orders the clusters:

1. The first `tmax` positions receive the clusters with the most user bins (partial sort).
2. The rest is sorted by the cardinality of each cluster's largest user bin, so that clusters with user bins of similar
   size are assigned after each other and small ones last.

If there are split bins, only `tmax - number_of_split_tbs` technical bins are left for merged bins and seeded. The
clusters at positions `[tmax - number_of_split_tbs, tmax)` are then assigned first, in order of their number of user
bins, not by cardinality as step 2 intends.

* **Variant A** (chopper): partial sort of the first `tmax` clusters.
* **Variant B**: partial sort of the first `number_of_remaining_tbs = tmax - number_of_split_tbs` clusters.

The variants differ only if there are split bins. Because the second calibration pass of `partition_user_bins` depends
on the first pass, a difference can remain even if the final partitioning has no split bins.

**Decision: keep variant A.** No systematic difference was found (see [Results](#results)).

### Method

`measure.cpp` generates one synthetic input, computes its sketches once and lays it out with both variants. The variant
is switched at run time by `chopper_variant_b`, which `variant_b.patch` adds to a copy of `partition_user_bins.cpp`.
With `chopper_variant_b = false`, the patched copy produces the same layouts as chopper.

**Inputs.** Each input imitates a collection of related genomes:

* User bins belong to families. Family k-mer pools have log-normal sizes (`median`, `sigma`), clamped to `[min, max]`.
* A user bin contains a random 60 % to 98 % of its family's pool, plus `(1 - share) * pool + 500` unique k-mers.
* `families = n / family_divisor`. `family_divisor = 1` means every user bin has its own family.
* `tmax` is chosen as by the command line interface (`seqan::hibf::config::validate_and_set_defaults`).

The inputs are deterministic. Running an input twice gives identical results.

| Grid | n | family_divisor | sigma | seeds | median, min, max (k-mers) | tmax | inputs |
|---|---|---|---|---|---|---|---|
| small | 400, 1000 | 1, 4, 25 | 0.7, 1.4 | 3 | 8000, 5000, 400000 | 64 | 36 |
| large | 10000, 25000, 50000 | 5, 50 | 1.4, 2.0 | 3, 2, 1 | 4000, 2500, 2000000 | 128, 192, 256 | 24 |

The minimum size is due to the MinHash sketches: computing them can fail for user bins with fewer than about 2000
k-mers.

**Metrics** (lower is better), per variant:

* `topmax`: the largest FPR-corrected technical bin of the top-level partitioning. It determines the size of the
  top-level IBF.
* `size`: expected total size of the HIBF in bytes (`hibf_statistics::total_hibf_size_in_byte`).
* `cost`: expected query cost of the HIBF (`hibf_statistics::expected_HIBF_query_cost`).

Also recorded: `s`, the number of top-level technical bins of user bins that span at least two technical bins, and
`same_top`, whether both variants produce the same top-level partitioning.

**Not covered.** A merged bin is laid out recursively with the fast layout only if it has at least `64 * tmax` user
bins. In the small grid this is impossible (`n < 64 * tmax`), so the lower levels were laid out by the DP algorithm of
the HIBF library and only the top-level partitioning differs between the variants. For the large grid, this was not
checked. All inputs are synthetic.

### Results

Ratios B/A over the inputs whose top-level partitionings differ. p-values are two-sided; the Wilcoxon signed-rank test
is on log(B/A). The per-input values are in `results_small.tsv` and `results_large.tsv`; `analyse.py` prints them.

**Small grid.** The top-level partitionings differ in 21 of 36 inputs.

| Metric | B better | B worse | equal | geo-mean B/A | range B/A | sign test p | Wilcoxon p |
|---|--:|--:|--:|--:|--:|--:|--:|
| largest top-level technical bin | 1 | 2 | 18 | 1.0004 | 0.979 - 1.031 | 1.00 | 0.75 |
| expected HIBF size | 7 | 10 | 4 | 0.9996 | 0.982 - 1.008 | 0.63 | 0.58 |
| expected query cost | 11 | 6 | 4 | 0.9988 | 0.971 - 1.015 | 0.33 | 0.52 |

**Large grid.** The top-level partitionings differ in only 9 of 24 inputs, and in none with 50000 user bins: with many
small user bins, the split threshold (joint size / `tmax`) is high and few user bins are split.

| Metric | B better | B worse | equal | geo-mean B/A | range B/A | sign test p | Wilcoxon p |
|---|--:|--:|--:|--:|--:|--:|--:|
| largest top-level technical bin | 4 | 2 | 3 | 0.9806 | 0.905 - 1.012 | 0.69 | 0.16 |
| expected HIBF size | 3 | 6 | 0 | 1.0062 | 0.997 - 1.034 | 0.51 | 0.25 |
| expected query cost | 3 | 5 | 1 | 1.0119 | 0.995 - 1.089 | 0.73 | 0.46 |

The largest effect is a single input (n = 10000, 200 families, sigma = 2.0, seed 3, s = 76): variant B shrinks the
largest top-level technical bin by 6.5 %, but the HIBF grows by 3.5 % and the query cost by 8.9 %. It dominates the
geometric means of the large grid.

Over both grids (30 inputs that differ), B is better in size for 10 inputs and worse for 16 (sign test p = 0.33), and
better in query cost for 14 and worse for 11 (p = 0.69). The differences go in both directions, are mostly below 1 %,
and no direction is significant. There is no measurable reason to change variant A.

### Reproduce

Requires CMake, a C++ compiler supported by chopper, `patch` and Python 3 (SciPy optional, for the Wilcoxon test).

```bash
CXX=clang++-23 ./build.sh /tmp/ppcv                  # configures a chopper Release build, builds /tmp/ppcv/measure
cd /tmp/ppcv && /path/to/this/directory/run.sh ./measure  # results_small.tsv, results_large.tsv (about 30 min)
/path/to/this/directory/analyse.py results_small.tsv results_large.tsv
```

`variant_b.patch` applies to `src/layout/partition_user_bins.cpp` as of the commit that added this directory. If the
file changes, the patch may need to be adapted.

`results_small.tsv` was produced by `run.sh`. `results_large.tsv` comes from an earlier run of the same input generator
before it was moved into `measure.cpp`; the first 12 of its 24 rows were checked against `run.sh` and are identical. The
timing columns (`t_sketch`, `t_A`, `t_B`, in seconds) depend on the machine.

### Environment

* chopper `e6fb149` (`final_fast_layout` branch), hibf `7f252fd`, Release build (`-O3 -DNDEBUG`).
* Compiler: Debian clang version 23.1.2, libstdc++.
* OS: Linux 7.2.8-WSL2-STABLE x86_64.
* CPU: AMD Ryzen 9 9950X, 16 cores, 32 threads. Memory: 49 GiB.
* Threads: 8 (small grid), 32 (large grid).
