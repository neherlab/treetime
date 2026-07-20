# Skyline optimization aborts on node-time inversions instead of degrading gracefully

`optimize_skyline` collects per-edge coalescent data via `collect_coalescent_edges` ([packages/treetime/src/coalescent/skyline.rs#L106](../../packages/treetime/src/coalescent/skyline.rs#L106)), which returns a hard error when any edge has `child_time < parent_time` ([packages/treetime/src/coalescent/edge_data.rs#L85-L90](../../packages/treetime/src/coalescent/edge_data.rs#L85-L90)). The error propagates through the final skyline re-optimization ([packages/treetime/src/timetree/pipeline.rs#L343-L361](../../packages/treetime/src/timetree/pipeline.rs#L343-L361)) and aborts the whole run.

The constant-Tc path handles the identical inversion gracefully: `optimize_tc` catches the error and falls back to the initial Tc with a warning ([packages/treetime/src/timetree/pipeline.rs#L234-L254](../../packages/treetime/src/timetree/pipeline.rs#L234-L254), [packages/treetime/src/timetree/pipeline.rs#L300-L318](../../packages/treetime/src/timetree/pipeline.rs#L300-L318)). Skyline does not, so `--coalescent-skyline` and `--coalescent-skyline --confidence` abort on inversion-prone datasets (e.g. flu/h3n2/200) while `--coalescent-opt` and `--coalescent=<v>` complete.

## Reproduction

```
timetree --tree=data/flu/h3n2/200/tree.nwk --dates=data/flu/h3n2/200/metadata.tsv \
  --aln=data/flu/h3n2/200/aln.fasta.xz --coalescent-skyline --n-skyline=5 --max-iter=3
```

Fails with `Failed to re-optimize skyline coalescent model` wrapping `Coalescent edge has child older than parent`.

## Notes

The inversions themselves stem from unconstrained marginal node-time selection, tracked in [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md), whose contract is undecided. This issue is narrower: the skyline path should treat an inversion the same way the constant-Tc path does (skip/degrade with a warning) rather than aborting, so the two coalescent code paths behave consistently. The graceful policy chosen for inversions must match the resolution of the marginal-node-time issue.

## Related issues

- [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md)
