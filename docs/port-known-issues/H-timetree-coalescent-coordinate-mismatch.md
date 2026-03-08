# Coalescent contributions use TBP coordinates, backward pass uses calendar time

[`compute_coalescent_contributions()`](../../packages/treetime/src/commands/timetree/coalescent/coalescent.rs#L61) (`#compute_coalescent_contributions`) produces distributions in time-before-present (TBP) coordinates, but the [backward pass](../../packages/treetime/src/commands/timetree/inference/backward_pass.rs#L45-L51) operates in calendar time. The two coordinate systems produce non-overlapping domains, so distribution multiplication returns `Distribution::Empty`. The coalescent prior has zero effect on node times.

## Root cause

v0 operates entirely in TBP coordinates: date constraints are converted to TBP
via `date2dist.get_time_before_present()` before entering the distribution
pipeline. The coalescent model, lineage counts, and merger rates all share the
same TBP domain. Distribution multiplication works because both operands share
the same axis.

v1 stores date constraints as calendar year fractions (e.g., 2005.5). The
backward pass propagates distributions in calendar time (e.g., domain
[1995, 2013]). The coalescent module converts calendar events to TBP before
computing contributions (e.g., domain [0, 18]). The two domains do not overlap.

Multiplication path in the backward pass:

1. Internal node `result` is initialized to the coalescent contribution
   (TBP domain, e.g., [0, 18])
2. First child message arrives in calendar time (e.g., [1995, 2013])
3. `distribution_multiplication` checks domain overlap:
   `overlap_min = max(0, 1995) = 1995`, `overlap_max = min(18, 2013) = 18`
4. Since `overlap_min >= overlap_max`, returns `Distribution::Empty`
5. All subsequent child message multiplications also return Empty
   (Empty absorbs in multiplication)

## Impact

- All internal node time distributions become `Distribution::Empty` when
  coalescent is enabled
- The forward pass multiplies parent messages with Empty subtree distributions,
  producing Empty results
- Internal nodes get `time = None` (Empty distributions have no `likely_time`)
- Leaf node times are unaffected (date constraints bypass the backward pass)
- The coalescent prior has no influence on divergence time inference

## Test gap

[`test_runner_coalescent_completes()`](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_runner_coalescent.rs#L37) (`#test_runner_coalescent_completes`) only checks that the pipeline does not panic and that extracted node times are finite. It uses `extract_node_times()` which filters out nodes with `time = None`, so internal nodes with Empty distributions are silently excluded. The test passes because leaf times are valid.

[`test_gm_coalescent()`](../../packages/treetime/src/commands/timetree/coalescent/__tests__/test_gm_coalescent.rs#L42) (`#test_gm_coalescent`) validates that contribution values match v0 snapshots in TBP space. This test is correct for what it tests (contribution computation), but does not cover integration with the backward pass.

## Fix direction

The coalescent contributions must be expressed in calendar time before entering
the backward pass. One approach: in `compute_coalescent_contributions`, shift
the DistributionFormula domain from TBP to calendar time by replacing `t` with
`present_time - t` in the evaluation closure and swapping `t_min`/`t_max` to
calendar coordinates.
