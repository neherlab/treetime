# Coalescent subsystem time notation conflict

## Summary

Signed time (negative = past) coexists with time-before-present (TBP, positive = past) in the same coalescent subsystem, creating ambiguity about the direction of time in integration bounds and distribution evaluations.

## Details

`packages/treetime/src/coalescent/`

Two conventions are used simultaneously:

1. Signed time (`coalescent.rs:34:`): past is negative t, present is 0, future is positive. Used for node dates and event ordering.

2. Time before present / TBP (`integration.rs:42-43:`): I(0) = 0 at the present, the integral increases going into the past. Positive values represent time depth.

The coalescent integral I(t) = integral from 0 to t of lambda(s) ds uses the TBP convention internally, but the integration bounds come from node times in the signed convention. The conversion between conventions is implicit (negation) and not documented at function boundaries.

## Required investigation

- Audit all functions in the coalescent module for which convention they expect as input and produce as output
- Verify that the implicit sign flip at the signed-to-TBP boundary is correct in all callers
- Determine if any callers pass signed times where TBP is expected (or vice versa)
- Compare against v0 Python implementation to verify convention consistency

## Impact

- If a sign flip is missed at any boundary, the coalescent integral is evaluated backward (present-to-past vs past-to-present)
- Wrong integration direction produces wrong coalescent costs, biasing time estimates
- The implicit convention makes the code fragile: any new caller must know the undocumented sign convention

## Related issues

- Source: [kb/issues/M-coalescent-time-notation-conflict.md](../issues/M-coalescent-time-notation-conflict.md) -- delete after full resolution
