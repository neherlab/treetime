# Timetree --keep-root rejects negative clock rate

When `--keep-root` is used, rerooting is disabled. If the input tree's root yields a negative clock rate from regression, v1 errors with "Estimated clock rate is non-positive" at [`packages/treetime/src/clock/clock_model.rs#L137`](../../packages/treetime/src/clock/clock_model.rs#L137). v0 proceeds with the negative rate and produces output.

## v0 behavior

v0 succeeds on flu/h3n2/20 with `--keep-root`, reporting rate -3.380e-03. It produces a timetree, molecular clock file, ancestral sequences, and auspice JSON. The negative rate is recorded in `molecular_clock.txt` and the tree is still time-scaled (though the interpretation is questionable).

## v1 behavior

v1 `ClockModel::from_regression()` at `clock_model.rs:136` checks `regression.clock_rate <= 0.0` and returns an error. This check fires during timetree pipeline before any inference begins.

## Affected datasets

- flu/h3n2/20 with `--keep-root`
- flu/h3n2/200 with `--keep-root`

Both are in the smoke test "Known failures" section (`dev/run-smoke-tests` lines 316-317).

## Analysis

The negative rate check is reasonable for the default rerooting path -- if all root positions yield negative rates, the data lacks temporal signal. But with `--keep-root`, the user explicitly opted to keep the input root. Rejecting a negative rate contradicts that intent.

v0's approach (proceed with negative rate) preserves user agency. The rate is reported in output, so the user can evaluate whether results are meaningful.

## Fix direction

When `--keep-root` is active, either:

- Skip the positive-rate check and proceed (v0 parity)
- Warn about negative rate and proceed
- Allow `--clock-rate` to override as the error message suggests, but also allow `--keep-root` to imply acceptance of the regression result

## Related issues

- [N-clock-regression-all-negative-rate.md](N-clock-regression-all-negative-rate.md) -- clock command also rejects all-negative rates, but that path reroots first
