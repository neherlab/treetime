# Guard non-positive Tc in coalescent contributions (yields NaN)

`fn compute_internal_contribution_single()` has no guard on non-positive Tc. When Tc(t) <= 0, lambda(t) goes negative or infinite, and `ln()` yields NaN or negative infinity. Callers do not validate Tc positivity before invoking.

v1: [`packages/treetime/src/coalescent/contributions.rs#L117-L178`](../../packages/treetime/src/coalescent/contributions.rs#L117-L178)

## Impact

Silent NaN/inf propagation in production builds for edge-case inputs. Degenerate inputs (non-positive Tc) not rejected at boundaries.

## Related issues

- Source: [N-numerical-stability-magic-constants.md](../issues/N-numerical-stability-magic-constants.md) -- delete after full resolution
