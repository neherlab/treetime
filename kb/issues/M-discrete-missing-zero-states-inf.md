# DiscreteNodeData::missing(0) produces inf-filled profile

`DiscreteNodeData::missing(n_states)` at [packages/treetime/src/representation/payload/discrete.rs#L24-L31](../../packages/treetime/src/representation/payload/discrete.rs#L24-L31) computes `1.0 / n_states as f64` to build a uniform prior profile. When `n_states == 0`, this evaluates to `1.0 / 0.0 = inf`, producing an `Array1` filled with `inf` values. No precondition check guards against zero states.

## Impact

Production code in the mugration command uses `n_states >= 1` (derived from trait alphabets that always have at least one state), so this path is not reachable in normal operation. The defect is a data-model soundness issue: the constructor accepts invalid input and produces a structurally unsound object rather than returning an error.

## Affected code

- Constructor: [packages/treetime/src/representation/payload/discrete.rs#L24-L31](../../packages/treetime/src/representation/payload/discrete.rs#L24-L31)
- Consumers: [packages/treetime/src/commands/mugration/discrete_marginal.rs](../../packages/treetime/src/commands/mugration/discrete_marginal.rs) (passes `partition.n_states()`, always >= 1)

## Fix

Add a precondition check: `if n_states == 0 { return make_internal_error!("n_states must be >= 1"); }` and change the return type to `Result<Self, Report>`. Callers derive `n_states` from alphabet size and never pass zero, so the error path is defensive.
