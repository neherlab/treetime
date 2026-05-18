# Guard zero states in discrete partition uniform_profile constructor

`uniform_profile(n_states)` at [packages/treetime/src/partition/marginal_discrete.rs#L165](../../packages/treetime/src/partition/marginal_discrete.rs#L165) computes `1.0 / n_states as f64` to build a uniform prior profile. When `n_states == 0`, this evaluates to `1.0 / 0.0 = inf`, producing an `Array2` filled with `inf` values. No precondition check guards against zero states.

## Impact

Production code in the mugration command validates `n_states < 2` at [packages/treetime/src/commands/mugration/run.rs#L186-L190](../../packages/treetime/src/commands/mugration/run.rs#L186-L190) and returns an error, so this path is not reachable in normal operation. The defect is a data-model soundness issue: the constructor accepts invalid input and produces a structurally unsound object rather than returning an error.

## Affected code

- Constructor: [packages/treetime/src/partition/marginal_discrete.rs#L165](../../packages/treetime/src/partition/marginal_discrete.rs#L165)
- Guard: [packages/treetime/src/commands/mugration/run.rs#L186-L190](../../packages/treetime/src/commands/mugration/run.rs#L186-L190)

## Fix

Add a precondition check in `PartitionMarginalDiscrete::new()` or in `uniform_profile()`: validate `n_states >= 1`. Low priority given the existing upstream guard.

## Related issues

- Source: [M-discrete-missing-zero-states-inf.md](../issues/M-discrete-missing-zero-states-inf.md)
