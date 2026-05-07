# Inference forward/backward pass asymmetries

## Summary

Four asymmetries between the backward and forward marginal passes: normalization strategy differs, graph mutation through interior mutability during read traversal, branch-length floor applied inside GTR inference, and clock model rebuilt unnecessarily.

## Instances

### Silent normalization asymmetry backward/forward pass

In the dense forward pass (`marginal_dense.rs:363:`), `msg_to_child` for each edge is computed as `node_data.profile.dis / safe_child` where `safe_child` is `msg_from_child.dis` clamped to `MIN_POSITIVE` (`marginal_dense.rs:424-428:`). The result is normalized via `normalize_inplace`.

The intermediate node profile at line 363 is built by multiplying `msg_to_parent.dis * msg_child` in probability space without per-step normalization across parent edges, unlike the backward pass which normalizes via `normalize_from_log` after combining child messages in log space.

This means the forward pass accumulates probability-space products (risk of underflow for deep trees) while the backward pass accumulates log-space sums (numerically stable).

### gather_clock_regression_results mutates graph via interior mutability

`packages/treetime/src/commands/clock/rtt.rs:42:`

Mutates node `div` (divergence) through a read-path graph traversal. The function signature suggests a read-only operation (gathering results) but has write side effects. This makes the function non-idempotent and its ordering relative to other graph operations significant.

### fix_branch_length clamp inside GTR inference

`packages/treetime/src/gtr/infer_gtr/dense.rs:192:`

Applies the marginal-pass branch-length floor during mutation counting for GTR inference. This means GTR parameter estimation sees clamped branch lengths rather than raw values, biasing the rate matrix toward shorter-branch statistics.

### run_refinement_iteration always rebuilds clock model

`packages/treetime/src/commands/timetree/refinement.rs:97:`

Rebuilds the clock model even when nothing moved (`n_diff==0 && n_resolved==0`). The clock model estimation involves regression over all dated tips, which is wasted computation when no dates changed. More importantly, floating-point non-determinism in the regression can introduce small perturbations that prevent clean convergence.

## Impact

- Forward pass underflow risk for deep trees with many siblings
- Hidden write side effects in nominally read-only functions
- GTR inference biased by branch-length clamping
- Unnecessary clock model rebuilds can prevent convergence detection

## Related issues

- Source: [M-inference-forward-backward-asymmetry.md](../issues/M-inference-forward-backward-asymmetry.md) -- delete after full resolution
