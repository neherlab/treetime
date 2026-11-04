# Smoke wall-time comparison mixes the load of concurrent cases

## Problem

`dev/smoke` runs as many cases at the same time as the container has CPUs (`fn default_jobs()`, the CPU quota of the container), and each case of a command with `--jobs` runs with the default thread count of TreeTime, so the cases compete for the same cores. The report then compares the wall time of each case with the wall time of the same case in the baseline snapshot ("Slowest cases", `+N% vs baseline`), although the two snapshots ran at different times, with different sets of cases running beside each case.

The wall time of a case therefore measures the load of its neighbors as much as the binary. A few long cases in one snapshot make the short cases of that snapshot slower, and the comparison reports this as a change of the binary.

## Evidence

In a smoke comparison of two binaries on a build host with 8 concurrent cases, the dense `optimize` cases of the candidate snapshot came out 1.5 to 4 times slower than in the baseline. Run alone, both binaries took the same time: `optimize` on `data/dengue/20` took 0.35 s with each binary, and on `data/ebola/100` 3.9 s and 3.2 s. During the candidate snapshot, three other cases held cores for about 120 s each.

## Fix direction

Make the timing comparison independent of the neighbors of a case, for example by one of these:

- Report the CPU time of each case (`result.json` already records it) instead of, or beside, the wall time
- Run the timed comparison with `--jobs 1`, or give each case a fixed thread count whose sum stays within the CPUs
- Mark wall-time changes as not comparable when the two snapshots ran different case sets at the same time

## Validation

Two smoke runs of the same binary against each other report no wall-time change above the run-to-run noise for the short cases, also when the selected case sets differ.
