# Timetree `--allow-negative-rate` has no v0 counterpart

v1 `timetree` accepts `--allow-negative-rate` and relaxes the positive-rate requirement of its root searches with it. v0 has the flag only on the `clock` command, and every reroot of the v0 timetree requires a positive rate. With the flag, a v1 timetree can therefore choose a different root than v0, and it reaches the outlier filter with a clock fit that v0 cannot produce. No decision approves the deviation.

## v0

- The flag is registered on the `clock` parser only ([packages/legacy/treetime/treetime/argument_parser.py#L447-L458](../../packages/legacy/treetime/treetime/argument_parser.py#L447-L458)) and read in `fn estimate_clock_model()` ([packages/legacy/treetime/treetime/wrappers.py#L984](../../packages/legacy/treetime/treetime/wrappers.py#L984)). The timetree parser has no such flag
- `fn TreeTime.run()` reroots through `fn TreeTime.reroot()` without `force_positive` ([packages/legacy/treetime/treetime/treetime.py#L260](../../packages/legacy/treetime/treetime/treetime.py#L260), [packages/legacy/treetime/treetime/treetime.py#L279](../../packages/legacy/treetime/treetime/treetime.py#L279)), and the default is `force_positive=True` ([packages/legacy/treetime/treetime/treetime.py#L542](../../packages/legacy/treetime/treetime/treetime.py#L542)). The clock filter of the timetree run reroots the same way ([packages/legacy/treetime/treetime/treetime.py#L486-L489](../../packages/legacy/treetime/treetime/treetime.py#L486-L489))

## v1

- CLI flag of timetree: `allow_negative_rate` in [packages/app-commands/src/commands/timetree/args.rs#L344-L346](../../packages/app-commands/src/commands/timetree/args.rs#L344-L346)
- Initial clock estimate: `fn estimate_initial_clock()` passes `!params.allow_negative_rate` as the positive-rate requirement ([packages/treetime/src/timetree/pipeline.rs#L285](../../packages/treetime/src/timetree/pipeline.rs#L285))
- Reroots of the pre-loop: `fn reroot()` passes the same value ([packages/treetime/src/timetree/pre_loop.rs#L136](../../packages/treetime/src/timetree/pre_loop.rs#L136)). The first reroot gives the clock fit that the outlier filter uses ([packages/treetime/src/timetree/pre_loop.rs#L43-L53](../../packages/treetime/src/timetree/pre_loop.rs#L43-L53)), so the flag also reaches the outlier filter of timetree. The pre-filter of the `clock` command requires a positive rate in every case ([packages/treetime/src/clock/pipeline.rs#L107](../../packages/treetime/src/clock/pipeline.rs#L107)), as v0 does
- The final positive-rate check of timetree stays in place: a non-positive rate at time inference fails the run ([kb/decisions/timetree-rejects-negative-clock-rate.md](../decisions/timetree-rejects-negative-clock-rate.md))

> [!IMPORTANT]
> **Decision required.** The timetree flag is a v1 deviation from v0 without approval. Options:
>
> - **Remove the flag from timetree** for v0 parity: every timetree reroot requires a positive rate, and the flag stays on `clock` only. Users who want a negative-rate regression use the `clock` command, as the decision file already advises
> - **Keep the flag** and record the deviation in a decision file: state the effect on the root searches and on the outlier filter of timetree, and why the relaxed root search serves timetree when the final check still rejects a non-positive rate
>
> Evidence: the v0 parser and calls above; the v1 calls above; the decision file, whose section "`--allow-negative-rate` in timetree" describes the v1 flag
