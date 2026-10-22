# Timetree reports one weak-signal failure with different errors

> [!IMPORTANT]
> **Investigation required.** Only one of the errors is recorded with its code path. The validation results that list the failing runs are not in the checkout. Produce them with the timetree suite of [validation/README.md](../../validation/README.md): build the binary with `./dev/docker/run just build-release treetime`, then in `validation/` run `pixi run python run.py --suite timetree --set core` and `pixi run python evaluate.py --suite timetree`, which writes `results/timetree.tsv` with the TreeTime error message of each failed run in its `error` column. From that table, collect the messages, the inputs that produce them, and the code path of each. The recorded reproductions of the `NaN` and undated-node errors (`data/dengue/500`, `data/dengue/1000`, `data/lassa/L/20`, and `keep-skyline-fixed/const/rep01` at $L = 10^4$ of the timetree suite) complete at commit `dd2eb7428`, so these two errors need new reproductions.

## Problem

On simulated data with weak temporal signal, the same case fails with different errors, depending on the code path. One of them is internal: `Cannot divide by empty distribution` [packages/treetime-distribution/src/distribution_ops/divide.rs#L30](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L30), raised by the division that forms the cavity of a node in the forward pass [packages/treetime/src/timetree/inference/forward_pass.rs#L195](../../packages/treetime/src/timetree/inference/forward_pass.rs#L195). The same division raises a second internal error for two point distributions [packages/treetime-distribution/src/distribution_ops/divide.rs#L31-L34](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L31-L34). None of the messages tells the user that the date signal is weak.

## Reproduction

The property test `test_prop_timetree_weak_signal_succeeds_or_rejects_the_input` [packages/treetime/src/timetree/__tests__/test_prop_timetree_weak_signal.rs](../../packages/treetime/src/timetree/__tests__/test_prop_timetree_weak_signal.rs) reproduces one weak-signal failure. The test is ignored until this issue is resolved; run it with `just test-rs --run-ignored only test_prop_timetree_weak_signal`. Its regression file `packages/treetime/proptest-regressions/timetree/__tests__/test_prop_timetree_weak_signal.txt` replays the case first.

- Case: the standard timetree cell of `fn timetree_cell()` (constant $T_c = 5$, 5 batches of 4 tips at $\tau = 0, 0.5, 1, 1.5, 2$, random GTR) with seed `1654128765928366562`, $\mu = 10^{-4}$ and $L = 221$; `realized` branch lengths, exact tip dates, an estimated clock rate, `coalescent: CoalescentMode::Constant` (as `--coalescent opt`), `keep_root: true`, `max_iter: 2`
- Result: `OperationError::InferenceFailed` with `The root-to-tip regression gives no positive clock rate, and the Gaussian timetree none either: Gaussian clock needs substitutions on the dated part of the tree`
- Path: the regression rate is not positive, so `fn ClockReroot.into_clock_fit_or()` [packages/treetime/src/clock/clock_regression.rs#L151-L162](../../packages/treetime/src/clock/clock_regression.rs#L151-L162) calls the fallback `fn gaussian_start_rate()` [packages/treetime/src/timetree/init.rs#L182-L192](../../packages/treetime/src/timetree/init.rs#L182-L192). The Gaussian tree fails with a plain report [packages/treetime/src/clock/gaussian.rs#L82-L84](../../packages/treetime/src/clock/gaussian.rs#L82-L84), and `fn OperationError::classify()` [packages/treetime/src/error.rs#L38-L46](../../packages/treetime/src/error.rs#L38-L46) maps every report without `InputError` or `NoPositiveRateRoot` to `InferenceFailed`. The fallback runs at the start clock [packages/treetime/src/timetree/pipeline.rs#L293-L294](../../packages/treetime/src/timetree/pipeline.rs#L293-L294) and in every round [packages/treetime/src/timetree/round.rs#L306-L307](../../packages/treetime/src/timetree/round.rs#L306-L307); the message does not tell which one failed
- Expected: the run succeeds, or returns `OperationError::InvalidInput` with a message that names the weak date signal, as `NoPositiveRateRoot` does for the regression alone. The message here names the cause, but the error kind reports an internal failure
- Frequency: 1 of 8 cases of the test's strategy ($L$ from 200 to 500, $\mu$ from $\{10^{-4}, 3 \cdot 10^{-4}\}$)

## Related errors

- `Gaussian clock needs substitutions on the dated part of the tree`: the fallback start rate on data without substitutions on the dated part of the tree, classified as an internal error (see Reproduction)

- `Cannot divide by empty distribution`: division by an empty distribution
- `Cannot divide point by point: operation not well-defined`: division of two point distributions
- a `NaN` time distribution, raised with and without round and node context
- an undated node, reported as a stale coalescent model by the point-time lineage counts, which phase 2 of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md) removed; the expected lineage counts count an undated node as surviving at all times

## Fix direction

Unknown until the paths are traced. A node whose date the rest of the tree gives no probability gets an empty posterior; the forward pass keeps the date the node was given and warns ([packages/treetime/src/timetree/inference/forward_pass.rs#L53-L63](../../packages/treetime/src/timetree/inference/forward_pass.rs#L53-L63)).

## Validation

Each weak-signal input of the validation suite ends with the same result on every code path: one error that names the cause, or a completed run.
