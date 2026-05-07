# Add convergence failure and warning reporting in GTR inference

## Description

Two locations where failures are silently suppressed in GTR inference and substitution composition.

### infer_gtr_impl silently proceeds after non-convergence

`packages/treetime/src/gtr/infer_gtr/common.rs:157-158:`

Returns `Ok(...)` with a `warn!` log when GTR inference does not converge. No convergence flag in the result struct. Callers cannot distinguish a converged model from a non-converged one without parsing log output.

Fix: add a convergence flag to the result struct so callers can programmatically detect and handle non-convergence.

### debug_assert_eq! stripped in release (compose_substitutions)

`packages/treetime/src/seq/mutation.rs:100:`

`debug_assert_eq!(ps.qry(), cs.reff(), ...)` is stripped in release builds. A broken substitution chain (where the query state of the parent substitution does not match the reference state of the child substitution) silently produces incorrect mutation annotations.

Fix: replace `debug_assert_eq!` with a proper error return or a runtime check that runs in all builds. Silent data corruption in release builds is not acceptable.

## Impact

- Non-converged GTR models used without caller awareness
- Silent data corruption in release builds from unchecked substitution chains

## Related issues

- Source: [N-error-suppression-unwrap-or-defaults.md](../issues/N-error-suppression-unwrap-or-defaults.md) -- delete after full resolution
