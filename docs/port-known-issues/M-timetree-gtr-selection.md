# GTR model selection not implemented

The timetree command unconditionally uses JC69 at
[`run.rs#L71`](../../packages/treetime/src/commands/timetree/run.rs#L71)
regardless of the `--model` CLI argument. The GTR model selection from
`args.gtr` is not wired through to the timetree inference pipeline. All
timetree runs use equal base frequencies and equal substitution rates, which is
incorrect for datasets with non-uniform composition. The ancestral command
implements GTR model selection correctly.
