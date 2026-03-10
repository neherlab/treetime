# Optimize sparse marginal reconstruction crashes with NaN

The `optimize` command with `--dense=false` (sparse marginal mode) crashes on larger
datasets during iterative branch length optimization.

## Repro

```bash
./dev/docker/run ./dev/dev r treetime -- optimize --dense=false \
  --tree=data/ebola/100/tree.nwk --aln=data/ebola/100/aln.fasta.xz \
  --outdir=tmp/optimize/sparse/ebola/100
```

Affected datasets: ebola/100, ebola/362, flu/h3n2/200, lassa/L/200

Working datasets: ebola/20, dengue/20, flu/h3n2/20

## Root Cause

Two bugs combine to cause the crash:

### Bug 1: NaN in probability vector (primary)

Location: `packages/treetime/src/representation/partition/marginal_helpers.rs:69`

```rust
if (*vec.max()? < (1.0 - EPS) * vec_norm) || !all_states_equal {
```

The `vec.max()` from `ndarray_stats::QuantileExt` fails with "Undefined ordering
between a tested pair of values" when `vec` contains NaN.

Cause: Probability vector multiplication produces underflow (very small values
multiplied together -> 0), then division by zero or log(0) operations introduce
NaN/Inf.

### Bug 2: Missing node re-insertion on error (secondary)

Location: `packages/treetime/src/representation/partition/marginal_passes.rs`

```rust
pub fn process_node_backward(...) {
  let mut seq_info = partition.nodes.remove(&node.key).unwrap();  // REMOVE
  // ... processing ...
  let result = combine_messages(...);
  result?  // EARLY RETURN on error - node never re-inserted!
  // ...
  partition.nodes.insert(node.key, seq_info);  // never reached on error
}
```

When `combine_messages` fails, the function returns early without re-inserting
the node. This corrupts the partition state, causing subsequent operations to
fail when trying to access the missing node.

## Error Message

```
The application panicked (crashed).
Message: called `Option::unwrap()` on a `None` value
Location: packages/treetime/src/representation/partition/marginal_passes.rs:25
```

Or with proper error propagation:

```
Undefined ordering between a tested pair of values.
Location: packages/treetime/src/representation/partition/marginal_helpers.rs:69
```

## Related Issues

- [N-ancestral-dense-normalize-log-nan](N-ancestral-dense-normalize-log-nan.md) -
  similar numerical stability issue in dense mode
