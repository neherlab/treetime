# Clock regression panics on a missing branch length

A tree whose edges have no branch lengths crashes the `clock` command:

```
treetime clock --tree=<tree without lengths> --metadata=data/zika/20/metadata.tsv --output-all=<dir>
Message:  Encountered an edge without a weight
Location: packages/treetime/src/clock/clock_regression.rs:379
```

`fn edge_divergence` calls `branch_length.expect("Encountered an edge without a weight")` ([packages/treetime/src/clock/clock_regression.rs#L368-L380](../../packages/treetime/src/clock/clock_regression.rs#L368-L380)). Other parts of the same command treat a missing length as 0: the divergence pass in `fn gather_clock_regression_results` ([packages/treetime/src/clock/rtt.rs#L30](../../packages/treetime/src/clock/rtt.rs#L30)) and `fn edge_branch_length` in the clock filter ([packages/treetime/src/clock/clock_filter.rs#L124](../../packages/treetime/src/clock/clock_filter.rs#L124)).

v0 replaces a missing branch length with 0 when it loads the tree ([packages/legacy/treetime/treetime/treeanc.py#L349](../../packages/legacy/treetime/treetime/treeanc.py#L349)), so its regression runs.

## Impact

- Input without branch lengths ends in a panic instead of either a v0-compatible result or an error message
- The regression and the reported root-to-tip divergences handle the same input in two ways

## Related

- [N-error-suppression-unwrap-or-defaults.md](N-error-suppression-unwrap-or-defaults.md): lists the default-to-0 sites
