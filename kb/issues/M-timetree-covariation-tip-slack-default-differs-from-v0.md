# Timetree covariation tip-slack default differs from v0

With `--covariation` and no `--tip-slack`, v1 timetree uses a tip slack of 10 (`tip_slack.unwrap_or(10.0)` in `fn build_covariation_clock_params`, [packages/treetime/src/timetree/params.rs#L72](../../packages/treetime/src/timetree/params.rs#L72)). The v1 `clock` command uses 3.

v0 has two defaults, and the CLI one applies to both commands:

- Python API: `ClockTree.__init__` sets `self.tip_slack = ttconf.OVER_DISPERSION`, which is 10 ([packages/legacy/treetime/treetime/clock_tree.py#L98](../../packages/legacy/treetime/treetime/clock_tree.py#L98), [packages/legacy/treetime/treetime/config.py#L8](../../packages/legacy/treetime/treetime/config.py#L8))
- CLI: `--tip-slack` defaults to 3 ([packages/legacy/treetime/treetime/argument_parser.py#L180-L185](../../packages/legacy/treetime/treetime/argument_parser.py#L180-L185)), and `run_timetree` assigns it to `myTree.tip_slack` ([packages/legacy/treetime/treetime/wrappers.py#L420](../../packages/legacy/treetime/treetime/wrappers.py#L420)), overriding the API default

So `treetime timetree --covariation` in v0 uses 3, and v1 uses 10.

## Impact

The tip term of the branch variance is `tip_slack²/L²`, so the v1 default gives terminal branches about 11 times more variance than v0 (100/9). Tips then carry less weight in the covariation-aware regression, which changes the fitted clock rate and the rerooting. Magnitude on real data was not measured.

## Related

- [M-clock-covariation-variance-diverges-from-v0.md](M-clock-covariation-variance-diverges-from-v0.md): the `clock` command's variance
- [kb/tickets/test-build-covariation-clock-params.md](../tickets/test-build-covariation-clock-params.md): unit tests of the same function
