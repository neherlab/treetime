# Scale invariance of partition rates is tested for one rate step only

With several alignment files, `treetime optimize` and `treetime timetree` estimate one relative rate per file (`--partition-rates`, [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md)). Multiplying every rate by a factor and dividing every branch length by it leaves the likelihood unchanged, so the estimated rates and branch lengths should not depend on the scale of the starting rates.

`test_optimize_rate_step_does_not_depend_on_the_scale_of_the_starting_rates` [packages/treetime/src/optimize/__tests__/test_pipeline_partitions.rs#L129-L149](../../packages/treetime/src/optimize/__tests__/test_pipeline_partitions.rs#L129-L149) checks one rate step, from starting rates 1, 2 and 0.5, at `1e-7` for the rates and `1e-8` for the branch lengths. The full optimization loop agrees only to about $2 \cdot 10^{-8}$ absolute in the branch lengths, $8 \cdot 10^{-6}$ relative on one short `flu/h3n2/20` branch (measured 2026-10-18), and no test covers the loop.

Two parts of branch-length optimization use fixed lengths that do not scale with the rates:

- **Floors**: the lower end of the Brent bracket is at least $10^{-12}$, and the indel floor is $0.01$ times the length of one mutation of the alignment [packages/treetime/src/optimize/method_brent.rs#L81-L86](../../packages/treetime/src/optimize/method_brent.rs#L81-L86), [packages/treetime/src/optimize/zero_boundary.rs#L32-L38](../../packages/treetime/src/optimize/zero_boundary.rs#L32-L38)
- **Tolerance**: the Brent solver stops at a relative tolerance of $1.49 \cdot 10^{-8}$ plus an absolute tolerance of $10^{-5}$ [packages/treetime/src/optimize/method_brent.rs#L13-L15](../../packages/treetime/src/optimize/method_brent.rs#L13-L15)

> [!IMPORTANT]
> **Investigation required.** Whether the loop-level differences come only from these fixed lengths, and whether they stay below the solver tolerance on every dataset. Run the loop from starting rates 1, 2 and 0.5 on the split datasets of the partition tests and record the largest absolute and relative differences of rates and branch lengths.

## Fix direction

Add a loop-level test of the scale invariance with a tolerance derived from the solver tolerance, once the investigation explains the differences.
