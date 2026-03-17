# Optimize sparse marginal reconstruction crashes with NaN

RESOLVED. The crash is fixed. The sparse `combine_messages()` (`#combine_messages`) now uses log-space arithmetic with logsumexp normalization, matching the dense mode approach. The unnecessary node removal from partition maps in backward and forward passes is eliminated.

Remaining: the optimize command still produces `-inf`/NaN log likelihood from iteration 3 onward on larger datasets (ebola/100, flu/h3n2/200). This affects both sparse and dense modes equally and is caused by degenerate distributions from the forward pass division operation, not by the `combine_messages` underflow that caused the crash. See [N-ancestral-dense-normalize-log-nan](N-ancestral-dense-normalize-log-nan.md).

## Original Repro

```bash
./dev/docker/run ./dev/dev r treetime -- optimize --dense=false \
  --tree=data/ebola/100/tree.nwk --aln=data/ebola/100/aln.fasta.xz \
  --outdir=tmp/optimize/sparse/ebola/100
```

## Fix

Two changes in [packages/treetime/src/representation/partition/marginal_helpers.rs](../../packages/treetime/src/representation/partition/marginal_helpers.rs):

- `combine_messages()` (`#combine_messages`): replaced direct probability vector multiplication (`vec *= &child_prob`) with log-space addition (`log_vec += child_prob.mapv(ln)`), followed by `logsumexp_normalize()` (`#logsumexp_normalize`). This prevents floating-point underflow when combining child messages at each tree node.

- Removed `ndarray_stats::QuantileExt::max()` call that failed on NaN. The variable-position check now uses a direct fold over the normalized distribution.

One change in [packages/treetime/src/representation/partition/marginal_passes.rs](../../packages/treetime/src/representation/partition/marginal_passes.rs):

- `process_node_backward()` (`#process_node_backward`) and `process_node_forward()` (`#process_node_forward`): replaced `partition.nodes.remove()` + `partition.nodes.insert()` with shared borrows and `get_mut()`. The node removal was a borrow-checker workaround that corrupted partition state on early return. Since `partition.nodes` and `partition.edges` are separate struct fields, Rust's borrow splitter handles concurrent access without removal.

## Related Issues

- [N-ancestral-dense-normalize-log-nan](N-ancestral-dense-normalize-log-nan.md) - similar numerical stability issue in dense mode, plus the shared forward-pass division problem that causes `-inf`/NaN on larger datasets
