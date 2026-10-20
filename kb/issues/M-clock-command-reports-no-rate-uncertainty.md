# The `clock` command reports no rate uncertainty

`treetime clock` reports the clock rate without a standard deviation. The Gaussian timetree that gives `rate_std` runs only in `timetree`, and the `clock` command has no sequence length to run it. The decision file states that the rate uncertainty comes from the Gaussian timetree and that v1 reports it in place of the regression covariance, without a restriction to `timetree`, so the decision text and the `clock` command disagree. v0 `treetime clock --covariation` prints the rate and the root date with one standard deviation each; v1 `clock` has no counterpart.

## Reproduction

```bash
./dev/docker/run just r treetime clock --tree=data/flu/h3n2/20/tree.nwk --metadata=data/flu/h3n2/20/metadata.tsv --output-all=tmp/clock/flu/h3n2/20
jq '.stats' tmp/clock/flu/h3n2/20/clock.clock-model.json
```

The statistics hold `chisq` and `r_val` only. `treetime timetree` on the same data with `--alignment` writes `stats.estimated.rate_std` (2.91e-4). Giving `clock` the alignment changes nothing, because `clock` does not read it.

## v0

- With `--covariation`, `DateConversion.__str__` prints the rate with the standard deviation from the covariance of the covariation-weighted regression ([packages/legacy/treetime/treetime/utils.py#L24-L32](../../packages/legacy/treetime/treetime/utils.py#L24-L32)), and the clock wrapper prints the root date with its standard deviation ([packages/legacy/treetime/treetime/wrappers.py#L1015-L1019](../../packages/legacy/treetime/treetime/wrappers.py#L1015-L1019)). Without `--covariation` it prints neither
- The covariation variances need the sequence length, from `--sequence-length` or the alignment, which the clock parser accepts ([packages/legacy/treetime/treetime/argument_parser.py#L134-L147](../../packages/legacy/treetime/treetime/argument_parser.py#L134-L147), [packages/legacy/treetime/treetime/argument_parser.py#L445](../../packages/legacy/treetime/treetime/argument_parser.py#L445))

## v1

- The clock model of a regression carries no standard deviation: `fn ClockModel::from_regression_unchecked()` sets `rate_std: None` ([packages/treetime/src/clock/clock_model.rs#L73-L83](../../packages/treetime/src/clock/clock_model.rs#L73-L83))
- Only `timetree` runs the Gaussian timetree and attaches its standard deviation to the clock model ([packages/treetime/src/timetree/pipeline.rs#L375-L379](../../packages/treetime/src/timetree/pipeline.rs#L375-L379))
- The Gaussian timetree needs the sequence length (`fn GaussianTree::new()` [packages/treetime/src/clock/gaussian.rs#L32-L41](../../packages/treetime/src/clock/gaussian.rs#L32-L41)). `clock` has none: `--sequence-length` was removed together with `--covariation`, and `clock` accepts `--alignment` without reading it ([packages/app-commands/src/commands/shared/alignment.rs#L33-L37](../../packages/app-commands/src/commands/shared/alignment.rs#L33-L37); `fn run_clock()` [packages/app-commands/src/commands/clock/run.rs#L24-L93](../../packages/app-commands/src/commands/clock/run.rs#L24-L93))
- Decision text without the restriction to `timetree`: [kb/decisions/clock-ordinary-least-squares-regression.md#L3-L5](../decisions/clock-ordinary-least-squares-regression.md#L3-L5) and the v0 divergence "Rate uncertainty" ([kb/decisions/clock-ordinary-least-squares-regression.md#L166-L167](../decisions/clock-ordinary-least-squares-regression.md#L166-L167)). The paragraph "Rate uncertainty" ([kb/decisions/clock-ordinary-least-squares-regression.md#L49-L57](../decisions/clock-ordinary-least-squares-regression.md#L49-L57)) ties the value to an estimated rate and to `--confidence`, which only `timetree` has

> [!IMPORTANT]
> **Decision required.** `clock` has no rate uncertainty, and the decision text does not say so. Options:
>
> - **Report the Gaussian uncertainty in `clock`**: give `clock` a sequence length (the length of `--alignment`, or a `--sequence-length` flag), run the Gaussian timetree at the regression rate, and write `rate_std` as `timetree` does; without a sequence length, write none. A root-date standard deviation as in v0 is a separate choice: the Gaussian root time is biased young against the regression root ([kb/proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md#L181](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md#L181)), so its width does not belong to the reported root date without a check
> - **Keep `clock` without uncertainty**: restrict the decision text to `timetree` and record the v0 divergence (`clock --covariation` printed the rate and the root date with one standard deviation)
>
> Evidence: the v0 and v1 locations above; `clock` and `timetree` runs on `data/flu/h3n2/20` with builds of `claude-clock-exact-root` (`e5eff6d07`), `claude-grid-trim` (`a3be8e9ff`) and the landed `rust`, all with the same result

## Validation

- First option: `clock` with a sequence length writes the same `rate_std` as the Gaussian timetree at the regression rate of the same tree and dates; `clock` without a sequence length writes none
- Second option: the decision text names `timetree` as the only command with a rate uncertainty, and the v0 divergence list of the decision file states the `clock` difference
