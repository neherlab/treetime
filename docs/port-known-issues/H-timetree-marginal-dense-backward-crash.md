# Marginal dense backward pass crash on ebola

The `PartitionMarginalDense` code path crashes in `process_node_backward` on
the ebola_20 dataset. The CLI default (sparse marginal) works; only the dense
path is affected.

## Location

Known skip in validation:
[`timetree_validation.rs#L233`](../../packages/treetime/examples/timetree_validation.rs#L233)

## Repro

```bash
./dev/docker/run ./dev/dev E timetree_validation
# Output: ebola_20/marginal_dense: SKIPPED (crashes in process_node_backward)
```

## Related issues

- [Dense backward pass produces NaN for all-zero probability rows](N-ancestral-dense-normalize-log-nan.md)
  shares the dense backward pass code path
