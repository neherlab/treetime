# Timetree inference in input mode collapses internal-node dates to Empty

When running timetree inference with `--branch-length-mode=input`, the inference pipeline multiplies point distributions for internal nodes. Repeated multiplication of narrow point distributions without stabilization collapses internal-node date distributions to `Empty`, losing the date information entirely.

## Mechanism

In input branch-length mode, the branch distributions are delta-like (point masses). The backward pass multiplies child messages (each a point distribution) with the branch distribution. Floating-point underflow in the product of narrow distributions produces an effectively zero distribution, which is classified as `Empty`.

Once a node's distribution becomes `Empty`, it is excluded from the forward pass and from date assignment, producing `None` dates for internal nodes that should have well-defined times.

## Impact

Internal node dates are missing from the nexus output. The timetree output contains leaf dates (from the input metadata) but no internal node dates. This is tracked separately in [M-timetree-internal-dates-missing-input-bl.md](M-timetree-internal-dates-missing-input-bl.md) as a symptom; the root cause is the distribution collapse described here.

## Fix

Use log-space arithmetic for distribution multiplication, or use a stabilized product that detects underflow and rescales before classifying the result as `Empty`.

## Related

- [M-timetree-internal-dates-missing-input-bl.md](M-timetree-internal-dates-missing-input-bl.md): output symptom (missing internal dates in nexus)
- [M-ancestral-marginal-probability-space.md](M-ancestral-marginal-probability-space.md): related probability-space issue in dense marginal
