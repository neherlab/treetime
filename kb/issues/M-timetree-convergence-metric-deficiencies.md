# Timetree convergence metric deficiencies

Remaining defects in the timetree convergence tracking system: the reported total log-likelihood is
not comparable across iterations, it plays no part in the convergence decision, and sequence-diff
counting undercounts across topology changes.

The original first defect — the convergence check keying on `n_diff`, which does not measure what
the loop moves — is resolved; see
[timetree-convergence-on-node-times.md](../decisions/timetree-convergence-on-node-times.md).
`has_converged` now tests node-time movement against `NODE_TIME_TOLERANCE_YEARS`, with `n_diff`
retained only as the fallback when no node is dated on both sides of a round.

## Details

### Likelihood plays no part in the convergence decision

[`metrics.rs`](../../packages/treetime/src/timetree/convergence/metrics.rs)

`ConvergenceMetrics` carries `log_lh_seq`, `log_lh_pos`, `log_lh_coal` and `log_lh_total`, and
`has_converged` uses none of them. A tree whose times have settled below tolerance while the
likelihood is still moving is declared converged. This cannot be fixed independently of the next
two items: the total is not currently a quantity a stop rule could be written against.

### Total log-likelihood is a sum of different terms across iterations

[`metrics.rs`](../../packages/treetime/src/timetree/convergence/metrics.rs)

`log_lh_total` is `[log_lh_seq, log_lh_pos, log_lh_coal].into_iter().flatten().reduce(...)`. When a
component is `Some` in one iteration and `None` in another the number of summed terms changes, so
deltas across iterations compare different objectives. `log_lh_coal` in particular is `None` or
`NaN` whenever `collect_coalescent_edges` fails, which happens silently — see
[M-coalescent-edge-collection-nan-bypass-and-unreachable-fallback.md](M-coalescent-edge-collection-nan-bypass-and-unreachable-fallback.md).

### The coalescent term is not the objective the times were inferred under

[`likelihood.rs`](../../packages/treetime/src/timetree/convergence/likelihood.rs) evaluates
`compute_coalescent_total_lh` against live lineage counts, while node times are inferred under a
prior built from counts frozen before the loop
([timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md)).
This is deliberate — the statistic should describe the tree you actually have — but it means
`log_lh_coal` is a diagnostic rather than a term of the maximized objective, and is not guaranteed
monotone. Any likelihood-based stop rule has to say which of the two it is testing.

### count_sequence_changes underreports on topology changes

[`sequence_changes.rs`](../../packages/treetime/src/timetree/convergence/sequence_changes.rs)

Compares per-partition sequence maps by zipping keys present in both snapshots. Nodes present in
only one — removed or created by polytomy resolution — are counted as `prev_only` / `curr_only` and
logged, but contribute no diffs. Lower severity than before, since `n_diff` is now only the
fallback criterion, but the count is still reported and written to the tracelog.

## Impact

- Convergence can be declared while the objective is still moving.
- `log_lh_total` deltas are not usable as a convergence signal without first fixing composition.
- The tracelog's `n_diff` column undercounts in rounds that changed the topology.

## Related

- [kb/tickets/timetree-convergence-metric-deficiencies.md](../tickets/timetree-convergence-metric-deficiencies.md)
- [N-timetree-convergence-tolerance-vs-branch-grid.md](N-timetree-convergence-tolerance-vs-branch-grid.md)
