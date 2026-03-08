# GTR model selection not implemented

The timetree command unconditionally uses JC69 at
[`run.rs#L71`](../../packages/treetime/src/commands/timetree/run.rs#L71)
regardless of the `--model` CLI argument. The GTR model selection from
`args.gtr` is not wired through to the timetree inference pipeline. All
timetree runs use equal base frequencies and equal substitution rates, which is
incorrect for datasets with non-uniform composition. The ancestral command
implements GTR model selection correctly.

Partition initialization at
[`initialization.rs#L87-L105`](../../packages/treetime/src/commands/timetree/initialization.rs#L87-L105)
always calls `jc69(JC69Params::default())`. The FIXME comment at
[`run.rs#L71`](../../packages/treetime/src/commands/timetree/run.rs#L71)
confirms this is known.

## Related issues

- [Dead CLI flags in timetree](N-timetree-dead-cli-flags.md) `--gtr-params`
  also dead
- [gtr.json missing for input-BL mode](N-timetree-gtr-json-missing-input-bl.md)
  no GTR output at all in input-BL mode
