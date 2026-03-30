# Grid search does not include t=0 as candidate branch length

## Severity

Negligible

## Description

The grid search fallback in `run_optimize_mixed()` evaluates branch lengths on `linspace(0.1 * one_mutation, 1.5 * branch_length + one_mutation, 100)`. This grid never includes $t = 0$ as a candidate.

When the Newton path executes (second derivative < 0 at the starting point), the clamping logic at `branch_length = 0.0` correctly converges to zero. The issue only manifests when the grid search branch executes (second derivative >= 0), which is uncommon for typical phylogenetic data.

## Affected code

- `run_optimize_mixed()` grid search at [packages/treetime/src/commands/optimize/optimize_unified.rs#L301](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L301)

## Impact

For edges where the true optimal branch length is zero and the grid search branch executes, the optimizer assigns a small positive branch length (`0.1 * one_mutation`) instead of zero. This is a bias toward positive branch lengths in the rare case where the Hessian is non-negative at $t = 0$.

The impact is more visible for non-unimodal models (K80, HKY85, TN93, GTR) after the model-aware zero-branch shortcut was added, since these models bypass the derivative-at-zero shortcut and always fall through to Newton/grid search.

## v0 behavior

v0 uses Brent's method (`scipy.optimize.minimize_scalar`) which searches the full interval and can find $t = 0$ as optimal.

## Fix

Include $t = 0$ in the grid search candidate set, or prepend 0.0 to the grid array.

## Dependencies

None.
