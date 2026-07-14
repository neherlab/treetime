# gtr.json missing for --branch-length-mode=input

`--branch-length-mode=input` output contains only `timetree.json`, `timetree.nexus`, `timetree.nwk`. No `gtr.json` is written because GTR model initialization is inside the `BranchLengthMode::Marginal` branch at [`run.rs#L105-L114`](../../packages/treetime/src/commands/timetree/run.rs#L105-L114).

v0 always writes `sequence_evolution_model.txt` regardless of branch length mode.

## Related issues

- Source: [kb/issues/N-timetree-gtr-json-missing-input-bl.md](../issues/N-timetree-gtr-json-missing-input-bl.md) -- delete after full resolution
- [Missing output files compared to v0](../issues/N-timetree-missing-output-files.md) -- umbrella issue
