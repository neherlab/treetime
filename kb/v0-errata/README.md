# v0 Errata

Defects, oversights, and undocumented simplifications in v0 (Python) that v1 (Rust) does not reproduce. One file per erratum.

Distinct from:

- [known issues](../issues/README.md) - v1 defects and missing features
- [decisions](../decisions/README.md) - deliberate v1 design choices that differ from correct v0 behavior

## When to add an erratum

Create an entry when v1 diverges from v0 because v0 is wrong. Check this directory before matching v0 behavior.

Evidence standard (at least two):

- The function accepts a parameter for the correct behavior but the caller omits it
- Adjacent code (same module, same author) handles the case correctly
- The docstring or comments describe the correct behavior but code does not implement it
- The behavior contradicts the cited scientific reference

## Summary

| Domain     | Erratum                                                                                       | v0 impact                                                  | v1 status |
| ---------- | --------------------------------------------------------------------------------------------- | ---------------------------------------------------------- | --------- |
| GTR        | [TN93 model ignores kappa parameters](tn93-alphabet-mismatch.md)                              | TN93 behaves as JC69                                       | Correct   |
| Confidence | [date_uncertainty_due_to_rate default interval typo](confidence-interval-default-typo.md)     | Dead code path (callers pass explicit args)                | Correct   |
| Coalescent | [total_LH uses fixed multiplicity=2 for all edges](coalescent-total-lh-fixed-multiplicity.md) | Tc optimization wrong for polytomies                       | Correct   |
| Optimize   | [optimize_tree_marginal uses signed convergence check](optimize-signed-convergence-check.md)  | Conflates convergence with failure (dead path in practice) | Improved  |
