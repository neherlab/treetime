# Clock regression rejects a missing branch length that v0 reads as 0

A tree whose edges have no branch lengths stops the `clock` command with an error:

```
treetime clock --tree=<tree without lengths> --metadata=data/zika/20/metadata.tsv --output-all=<dir>
Error: Edge 2 has no branch length
```

`fn edge_divergence` returns this error when an edge has neither a branch length nor a clock length ([packages/treetime/src/clock/clock_regression.rs#L445-L458](../../packages/treetime/src/clock/clock_regression.rs#L445-L458)), and the reroot step of the clock tree does the same for an inverted edge ([packages/treetime/src/clock/reroot.rs#L341](../../packages/treetime/src/clock/reroot.rs#L341)). Other parts of the same command treat a missing length as 0 through `branch_length_or_zero()`: the root-to-tip divergences that the clock output reads ([packages/treetime/src/clock/pipeline.rs#L68](../../packages/treetime/src/clock/pipeline.rs#L68)) and the clock filter ([packages/treetime/src/clock/clock_filter.rs#L37](../../packages/treetime/src/clock/clock_filter.rs#L37)).

v0 replaces a missing branch length with 0 when it loads the tree ([packages/legacy/treetime/treetime/treeanc.py#L349](../../packages/legacy/treetime/treetime/treeanc.py#L349)), so its regression runs.

## Impact

- Input without branch lengths ends in an error instead of the v0 result
- The regression and the reported root-to-tip divergences handle the same input in two ways

> [!IMPORTANT]
> **Decision required.** The regression either reads a missing branch length as 0, as v0 and the divergence pass do, or keeps the error and the divergence pass and clock filter report the same error. The v0 behavior is the default reference; the error is a deviation until approved.

## Related

- [N-error-suppression-unwrap-or-defaults.md](N-error-suppression-unwrap-or-defaults.md): lists the default-to-0 sites
