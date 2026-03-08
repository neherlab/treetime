# Internal node dates missing at scale

Date annotations in nexus output decrease with tree size. flu/h3n2/20 has 37/37
annotations (100%), flu/h3n2/500 has 476/917 (52%). The positional likelihood
is `NaN` because `compute_positional_likelihood()` skips edges where parent or
child time is `None`.

## Root cause

The backward pass produces time distributions for internal nodes, but the
forward pass at
[`forward_pass.rs#L46`](../../packages/treetime/src/commands/timetree/inference/forward_pass.rs#L46)
calls `dist.likely_time()` which returns `None` for some distributions. These
nodes get `time = None` and no nexus annotation.

The problem scales with tree size: distribution computation becomes degenerate
(too narrow, too wide, or numerically unstable) for longer path lengths in
larger trees.

## Positional likelihood

`compute_positional_likelihood()` at
[`likelihood.rs#L31`](../../packages/treetime/src/commands/timetree/convergence/likelihood.rs#L31)
skips edges where parent or child time is `None`. When most internal nodes have
`None` time, the metric returns `None` (displayed as NaN).

## Compound with marginal CIs

When combined with `--time-marginal=only-final` on flu/h3n2/500, the
`confidence_intervals.tsv` file contains 476 entries (all tips), with zero
internal nodes. All tip CIs have zero width. CI estimation is useless at this
scale because internal node dates are dropped before marginal reconstruction
runs.

## Related issues

- [Clock rate frozen across iterations](H-timetree-clock-rate-frozen.md) makes
  the rate less accurate at scale, worsening distribution degeneration
- [Internal node dates missing for coalescent modes](M-timetree-internal-dates-missing-coalescent.md)
  similar symptom from a different cause (coalescent distributions)
- [Internal dates missing with bad fixed rate](M-timetree-internal-dates-bad-fixed-rate.md)
  similar mechanism triggered by mismatched clock rate
