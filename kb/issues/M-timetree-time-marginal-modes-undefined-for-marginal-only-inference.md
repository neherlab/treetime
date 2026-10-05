# `--time-marginal` modes are undefined for marginal-only time inference

v1 infers node times only by marginal inference ([kb/decisions/timetree-marginal-only-time-inference.md](../decisions/timetree-marginal-only-time-inference.md)). The `--time-marginal` values come from v0, where they choose between joint and marginal inference. In v1 they choose something else, and no entry defines what.

## v0 meaning

- `never` (v0 default `false`): joint time inference in every round
- `always` (`true`): marginal time inference in every round
- `only-final` and `confidence-only`: joint inference in every round, then one final marginal round that re-estimates the clock and overwrites the dates with the marginal peaks

Sources: [packages/legacy/treetime/treetime/argument_parser.py#L258-L265](../../packages/legacy/treetime/treetime/argument_parser.py#L258-L265), [packages/legacy/treetime/treetime/treetime.py#L218-L222](../../packages/legacy/treetime/treetime/treetime.py#L218-L222), [packages/legacy/treetime/treetime/treetime.py#L390-L396](../../packages/legacy/treetime/treetime/treetime.py#L390-L396)

## v1 behavior

The CLI accepts `never`, `always`, and `only-final` (`enum TimeMarginalModeCli` in [packages/app-commands/src/commands/timetree/args.rs](../../packages/app-commands/src/commands/timetree/args.rs)). Every round runs the same marginal inference for all three:

- `never` and `always` give identical dates
- `always` and `only-final` write `confidence_intervals.tsv` from the marginal posteriors; `never` does not ([kb/decisions/timetree-ci-output-ungated.md](../decisions/timetree-ci-output-ungated.md), `fn gather_results` in [packages/treetime/src/timetree/pipeline.rs](../../packages/treetime/src/timetree/pipeline.rs))
- `only-final` also runs one extra inference round after the loop, commits the clock branch lengths without damping, and runs a marginal update of the sequence partition; v0's final round has no ancestral update
- `--confidence` changes `never` to `only-final` when `--clock-std-dev` or `--covariation` is given, and warns otherwise (`fn compute_effective_time_marginal` in [packages/treetime/src/timetree/params.rs](../../packages/treetime/src/timetree/params.rs))

The flag name and the values `never` and `always` suggest a choice of inference method that no longer exists.

> [!IMPORTANT]
> **Decision required.** What should the timetree command offer instead of the v0 meaning of `--time-marginal`?
>
> - Keep the three values and define them by what they do in v1: `never` gives dates without interval output, `always` adds intervals from the regular rounds, `only-final` adds the final round. No code change; the flag name and the values `never` and `always` stay misleading
> - Replace the flag with options named after what they control (interval output, extra final round), and remove `--time-marginal`. This follows the project rule to remove superseded options and changes the CLI, the JSON schemas, and the analysis form
>
> Evidence: the v1 behavior above. The choice also decides whether `only-final` keeps its extra round, which no v1 entry justifies on its own.

## Validation

After the decision: CLI help text, `kb/features/timetree.md` section "Time Marginal Modes", and [kb/decisions/timetree-ci-output-ungated.md](../decisions/timetree-ci-output-ungated.md) describe the same behavior; a test per mode checks the dates and which outputs are written.
