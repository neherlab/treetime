# gtr.json missing for --branch-length-mode=input

`--branch-length-mode=input` output contains only `timetree.json`, `timetree.nexus`, `timetree.nwk`. No `gtr.json` is written because GTR model initialization is inside the `BranchLengthMode::Marginal` branch; the `BranchLengthMode::Input` branch returns `gtr: None` at [`pipeline.rs#L332-L338`](../../packages/treetime/src/timetree/pipeline.rs#L332-L338).

v0 always writes `sequence_evolution_model.txt` regardless of branch length mode.
