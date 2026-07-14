# Timetree convergence metric deficiencies

## Summary

Three defects in the timetree convergence tracking system: the convergence check ignores likelihood stagnation, topology changes cause undercounting of sequence diffs, and the total likelihood composition changes between iterations.

## Details

### has_converged ignores likelihood stagnation

`packages/treetime/src/timetree/convergence/metrics.rs:140-142:`

`fn has_converged()` returns true when `n_diff==0 && n_resolved==0` (no node date changes, no newly resolved dates). It does not check whether `lh_total` has plateaued. The `ConvergenceMetrics` struct carries `lh_seq`, `lh_pos`, `lh_coal`, and `lh_total` fields, but none are used in the convergence decision. A tree where dates are stable but likelihood is still improving (or worsening) is declared converged.

### count_sequence_changes underreports on topology changes

`packages/treetime/src/timetree/convergence/sequence_changes.rs:25-44:`

Compares per-partition sequence maps between iterations by zipping keys present in both `previous` and `current`. Nodes that exist only in `previous` (removed by topology change) or only in `current` (newly created) are counted as `prev_only` / `curr_only` but not as sequence diffs. The total diff count excludes all changes at nodes that were replaced during polytomy resolution or pruning.

### total_lh_reduce meaning changes with Some/None patterns

`packages/treetime/src/timetree/convergence/metrics.rs:67:`

`lh_total` is computed as `[lh_seq, lh_pos, lh_coal].into_iter().flatten().reduce(|acc, v| acc + v)`. When the coalescent prior is enabled in one iteration but not another (first iteration runs without coalescent, second adds it), the number of `Some` components changes. Comparing `lh_total` across iterations compares sums of different terms, making likelihood deltas meaningless as convergence indicators.

## Impact

- Premature convergence declaration when dates stabilize but likelihood has not plateaued
- After polytomy resolution, convergence check undercounts actual sequence changes
- Likelihood-based convergence criteria are unreliable when coalescent is enabled mid-run

## Related issues

- Source: [kb/issues/M-timetree-convergence-metric-deficiencies.md](../issues/M-timetree-convergence-metric-deficiencies.md) -- delete after full resolution
