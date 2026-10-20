# gtr.json missing for --branch-length-mode=input

With `--branch-length-mode=input`, the timetree pipeline builds no partitions and fits no GTR model. `fn initialize_branch_model()` returns an empty partition list and `BranchModel::Input` for this mode, and builds the partitions from the alignment only in marginal mode ([packages/treetime/src/timetree/pipeline.rs#L305-L316](../../packages/treetime/src/timetree/pipeline.rs#L305-L316)). The command therefore has no model to write to `timetree.gtr.json`.

The skip reason is misleading. `fn write_gtr_output()` receives the reason "no GTR model was fitted (provide an alignment with --alignment)" ([packages/app-commands/src/commands/timetree/run.rs#L217-L220](../../packages/app-commands/src/commands/timetree/run.rs#L217-L220)). With `--output-gtr`, the run fails with this message; with `--output-all`, the file is skipped with a debug message ([packages/app-output/src/output_plan.rs#L528-L541](../../packages/app-output/src/output_plan.rs#L528-L541)). In input mode the message tells the user to provide an alignment also when `--alignment` was given.

v0 writes `sequence_evolution_model.txt` whenever `--gtr infer` is set (the default), independent of the branch length mode ([packages/legacy/treetime/treetime/wrappers.py#L418](../../packages/legacy/treetime/treetime/wrappers.py#L418), [packages/legacy/treetime/treetime/wrappers.py#L511-L516](../../packages/legacy/treetime/treetime/wrappers.py#L511-L516)).

Every run with `--branch-length-mode=input` currently aborts before it writes outputs ([H-timetree-input-branch-lengths-abort-on-point-division.md](H-timetree-input-branch-lengths-abort-on-point-division.md)), so this gap becomes visible only after that fix.
