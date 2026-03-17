# gtr.json missing for --branch-length-mode=input

`--branch-length-mode=input` output contains only `timetree.json`,
`timetree.nexus`, `timetree.nwk`. No `gtr.json` is written because GTR model
initialization is inside the `BranchLengthMode::Marginal` branch at
[`run.rs#L69-L74`](../../packages/treetime/src/commands/timetree/run.rs#L69-L74).

v0 always writes `sequence_evolution_model.txt` regardless of branch length
mode.

## Related issues

- [Missing output files compared to v0](N-timetree-missing-output-files.md)
  umbrella issue
