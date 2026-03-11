# Iterative GTR inference not implemented for mugration

The mugration command creates a uniform GTR model at [`run.rs#L124-L129`](../../packages/treetime/src/commands/mugration/run.rs#L124-L129) with `W: None` (equal transition rates) and either uniform or weight-based equilibrium frequencies. The model is used directly for a single forward-backward reconstruction pass without iterative refinement.

v0 `reconstruct_discrete_traits()` (`#reconstruct_discrete_traits`) at [`wrappers.py#L653-L811`](../../packages/legacy/treetime/treetime/wrappers.py#L653-L811) performs iterative GTR inference:

1. Initial `infer_ancestral_sequences(infer_gtr=True)` estimates the rate matrix from observed transitions ([`wrappers.py#L785-L794`](../../packages/legacy/treetime/treetime/wrappers.py#L785-L794))
2. `optimize_gtr_rate()` optimizes the overall rate scalar via Brent minimization ([`wrappers.py#L795`](../../packages/legacy/treetime/treetime/wrappers.py#L795))
3. Five iterations of `infer_gtr()` + `optimize_gtr_rate()` re-estimate equilibrium frequencies and rates from the data ([`wrappers.py#L800-L802`](../../packages/legacy/treetime/treetime/wrappers.py#L800-L802))
4. Final `infer_ancestral_sequences(infer_gtr=False)` reconstructs with the refined model ([`wrappers.py#L807-L809`](../../packages/legacy/treetime/treetime/wrappers.py#L807-L809))

The iterative refinement shifts equilibrium frequencies away from uniform, changing posterior profiles at ambiguous internal nodes. Golden master tests confirm that 2/6 datasets (zika, lassa) produce identical trait assignments despite the difference, while 4/6 datasets (dengue, tb, rsv, mpox) diverge at internal nodes where the phylogeographic signal is weak.

## Related issues

- [GTR model selection not implemented](M-timetree-gtr-selection.md) same gap for timetree command
- [Iterative GTR for discrete traits](../port-algo-inventory/unimplemented.md#iterative-gtr-for-discrete-traits) algorithm description in unimplemented inventory
