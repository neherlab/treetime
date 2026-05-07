# Coalescent integration test midpoint-rule discretization error

`fn test_compute_integral_merger_rate_varying_tc_many_segments` at [packages/treetime/src/commands/timetree/coalescent/\_\_tests\_\_/test_integration.rs#L157](../../packages/treetime/src/commands/timetree/coalescent/__tests__/test_integration.rs#L157) computes the numerical integral of `1/(0.01 + 0.004t)` over [0, 10] using midpoint-rule segments and compares against the analytical result `250 * ln(5) = 402.36`.

Tolerance set to `1e-6`, test `#[ignore]`d. With 1000 segments the midpoint-rule error is 1.6e-4. O(h^2) convergence requires ~12500 segments for 1e-6. The production code is correct; the issue is the test's numerical quadrature accuracy.

## Fix

Replace midpoint-rule quadrature with a higher-order method (e.g. Simpson's rule, Gauss-Legendre) or use enough segments for O(h^2) convergence to reach 1e-6.
