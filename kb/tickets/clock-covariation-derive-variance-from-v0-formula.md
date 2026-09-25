# Derive the clock command's covariation variance from the v0 formula

Replace the app-layer variance construction of the `clock` command with the core function timetree already uses, so both commands compute the v0 covariation variance in one place.

## Current state

[packages/app-commands/src/commands/clock/run.rs#L57-L68](../../packages/app-commands/src/commands/clock/run.rs#L57-L68) builds `ClockVarianceParams` with `variance_factor = 2.0 / seq_len` and takes `seq_len` only from `--sequence-length`. `fn build_covariation_clock_params` in [packages/treetime/src/timetree/params.rs#L54-L85](../../packages/treetime/src/timetree/params.rs#L54-L85) implements the v0 formula (`1/L`, leaf offset `tip_slack²/L²`) and derives `L` from the alignment when one is given.

## Task

- Make `build_covariation_clock_params` reachable from `app-commands` (it is `pub(crate)` in `treetime`) and call it from the clock command, passing the clock command's alignment (when given), `--sequence-length` and `--tip-slack`
- Delete the local `overdispersion = 2.0` construction
- Keep the clock command's `--tip-slack` default at 3 (v0 CLI default, [packages/legacy/treetime/treetime/argument_parser.py#L180-L185](../../packages/legacy/treetime/treetime/argument_parser.py#L180-L185)). The timetree default is tracked separately in [kb/issues/M-timetree-covariation-tip-slack-default-differs-from-v0.md](../issues/M-timetree-covariation-tip-slack-default-differs-from-v0.md); if that issue is resolved first, one shared default serves both commands
- Warn when `--covariation` overrides explicitly given `--variance-factor`, `--variance-offset` or `--variance-offset-leaf`
- Test: `clock --covariation --alignment=<aln>` without `--sequence-length` succeeds and uses the alignment length; the variance parameters equal the v0 formula for a hand-computed $L$ and tip slack (oracle: [packages/legacy/treetime/treetime/clock_tree.py#L277-L285](../../packages/legacy/treetime/treetime/clock_tree.py#L277-L285))
- Run smoke on the clock cases with `--covariation`; changed clock rates in those cases are expected and must be explained by the factor change

## Related issues

- Source: [kb/issues/M-clock-covariation-variance-diverges-from-v0.md](../issues/M-clock-covariation-variance-diverges-from-v0.md) -- delete after full resolution
