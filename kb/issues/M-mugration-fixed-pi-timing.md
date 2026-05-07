# Mugration fixed_pi capture happens before pseudo-count smoothing

In `execute_mugration` at [packages/treetime/src/commands/mugration/run.rs#L216-L221](../../packages/treetime/src/commands/mugration/run.rs#L216-L221), `fixed_pi` is captured as `weights.as_ref().map(|_| pi.clone())` before `apply_pseudo_counts` is called. The captured `fixed_pi` preserves the raw weighted pi (pre-pseudo-count), matching v0's `fixed_pi=weights` intent.

The default pseudo-count `pc.unwrap_or(1.0)` at [run.rs#L251](../../packages/treetime/src/commands/mugration/run.rs#L251) is forwarded to `refine_gtr_iterative` regardless of whether the user set `--pc`. The CLI `--pc` arg at [args.rs#L40](../../packages/treetime/src/commands/mugration/args.rs#L40) has no `default_value`, leaving `pc: Option<f64>`. Users passing no `--pc` flag get `pc=1.0` inside refinement but `pc=None` (no smoothing) in `apply_pseudo_counts`, creating an inconsistency between the two pi-smoothing paths.

## Affected code

- `fixed_pi` capture: [packages/treetime/src/commands/mugration/run.rs#L216-L221](../../packages/treetime/src/commands/mugration/run.rs#L216-L221)
- Default pseudo-count: [packages/treetime/src/commands/mugration/run.rs#L251](../../packages/treetime/src/commands/mugration/run.rs#L251)
- CLI arg: [packages/treetime/src/commands/mugration/args.rs#L40](../../packages/treetime/src/commands/mugration/args.rs#L40)

## Fix

Unify the pseudo-count default: either use `default_value_t = 1.0` on the CLI arg and drop the `unwrap_or`, or pass `None` through both paths so that the same smoothing (or lack thereof) applies to both `apply_pseudo_counts` and `refine_gtr_iterative`.

## v0 comparison

v0 uses `pc=1.0` inside `infer_gtr` by default. The v1 behavior matches v0 for the refinement path but diverges for the initial `apply_pseudo_counts` path when no `--pc` is specified. Document this as an intentional change if the split is desired.
