# Timetree `--allow-negative-rate` has no v0 counterpart

v1 `timetree` accepts `--allow-negative-rate` and relaxes the positive-rate requirement of its root searches with it. v0 has the flag only on the `clock` command, and every reroot of the v0 timetree requires a positive rate. With the flag, a v1 timetree can therefore choose a different root than v0, and it reaches the outlier filter with a clock fit that v0 cannot produce. No decision approves the deviation.

## v0

- The flag is registered on the `clock` parser only ([packages/legacy/treetime/treetime/argument_parser.py#L447-L458](../../packages/legacy/treetime/treetime/argument_parser.py#L447-L458)) and read in `fn estimate_clock_model()` ([packages/legacy/treetime/treetime/wrappers.py#L984](../../packages/legacy/treetime/treetime/wrappers.py#L984)). The timetree parser has no such flag
- `fn TreeTime.run()` reroots through `fn TreeTime.reroot()` without `force_positive` ([packages/legacy/treetime/treetime/treetime.py#L260](../../packages/legacy/treetime/treetime/treetime.py#L260), [packages/legacy/treetime/treetime/treetime.py#L279](../../packages/legacy/treetime/treetime/treetime.py#L279)), and the default is `force_positive=True` ([packages/legacy/treetime/treetime/treetime.py#L542](../../packages/legacy/treetime/treetime/treetime.py#L542)). The clock filter of the timetree run reroots the same way ([packages/legacy/treetime/treetime/treetime.py#L486-L489](../../packages/legacy/treetime/treetime/treetime.py#L486-L489))

## v1

- CLI flag of timetree: `allow_negative_rate` of the shared clock filter arguments in [packages/app-commands/src/commands/shared/clock_filter.rs#L40-L48](../../packages/app-commands/src/commands/shared/clock_filter.rs#L40-L48)
- Initial clock estimate: `fn estimate_initial_clock()` passes `!params.allow_negative_rate` as the positive-rate requirement ([packages/treetime/src/timetree/pipeline.rs#L406](../../packages/treetime/src/timetree/pipeline.rs#L406))
- Reroots of the pre-loop: `fn reroot()` passes the same value ([packages/treetime/src/timetree/pre_loop.rs#L248](../../packages/treetime/src/timetree/pre_loop.rs#L248)). The first reroot gives the clock fit that the outlier filter uses ([packages/treetime/src/timetree/pre_loop.rs#L44-L50](../../packages/treetime/src/timetree/pre_loop.rs#L44-L50)), so the flag also reaches the outlier filter of timetree. The pre-filter of the `clock` command requires a positive rate in every case ([packages/treetime/src/clock/pipeline.rs#L109](../../packages/treetime/src/clock/pipeline.rs#L109)), as v0 does
- With a `least-squares` root, the pre-loop refits keep the positive-rate check: a non-positive rate there fails the run. The other clock fits take the rate of the Gaussian timetree in place of a non-positive one ([kb/decisions/timetree-runs-no-time-inference-on-a-non-positive-clock-rate.md](../decisions/timetree-runs-no-time-inference-on-a-non-positive-clock-rate.md))

> [!IMPORTANT]
> **Decision required.** The timetree flag is a v1 deviation from v0 without approval. Options:
>
> - **Remove the flag from timetree** for v0 parity: every timetree reroot requires a positive rate, and the flag stays on `clock` only. Users who want a negative-rate regression use the `clock` command, as the decision file already advises
> - **Keep the flag** and record the deviation in a decision file: state the effect on the root searches and on the outlier filter of timetree, and why the relaxed root search serves timetree when the final check still rejects a non-positive rate
>
> Evidence: the v0 parser and calls above; the v1 calls above; the decision file, whose section "`--allow-negative-rate` in timetree" describes the v1 flag

## Proposed implementation plan

The work waits for the decision above. One plan per option:

**Option: remove the flag from timetree**

1. Split `allow_negative_rate` out of the shared clock filter arguments in `packages/app-commands/src/commands/shared/clock_filter.rs#L40-L48` into the `clock` command arguments only, so `timetree --help` no longer lists it.
2. Remove `allow_negative_rate` from the timetree parameters and pass `true` as the positive-rate requirement in `fn estimate_initial_clock()` (`packages/treetime/src/timetree/pipeline.rs#L406`) and `fn reroot()` (`packages/treetime/src/timetree/pre_loop.rs#L248`); delete code that only served the relaxed path.
3. Remove the section "`--allow-negative-rate` in timetree" from `kb/decisions/timetree-runs-no-time-inference-on-a-non-positive-clock-rate.md` (needs user consent, as it is a decision file), and add the spelling to `packages/app-cli/src/cli/v0_args.rs` only if v0 timetree accepted it (it did not, so no entry).
4. Tests: a CLI parse test that `treetime timetree --allow-negative-rate` fails with clap's unknown-argument error (oracle: v0 timetree parser, `packages/legacy/treetime/treetime/argument_parser.py`); a timetree test on a dataset where the best root has a negative rate (the `data/dengue/100` scenario of `test_clock_dengue100.rs`) that asserts the same root as v0 (oracle: v0 `treetime` run on that dataset).
5. Validate with `just check`, `just test-rs`, `just generated-check` (CLI arguments and schemas change), and `just smoke`.

**Option: keep the flag**

1. With user approval, write the deviation into a decision file: the effect on the initial clock estimate, on the pre-loop reroots, and on the outlier filter of timetree (`packages/treetime/src/timetree/pre_loop.rs#L44-L50`), and the reason.
2. Update the `--help` text of the flag in `clock_filter.rs` to name the timetree effect on root search and outlier filter.
3. Tests: a timetree test on the dengue/100 scenario that asserts that the flag changes the chosen root and that the final check still rejects a non-positive rate (oracle: the `clock` command with `--allow-negative-rate` on the same tree, which chooses the same root).
4. Validate with `just check`, `just test-rs`, and `just generated-check`.
