# Fitch recurrence is not minimum parsimony on multifurcations

The sparse Fitch backward pass implements Fitch's <a id="gloss-use-1"></a>minimum-change method <sup>[1](#gloss-1)</sup> by intersecting all child state sets and falling back to their union when the common intersection is empty [packages/treetime/src/ancestral/fitch_sub.rs#L38-L72](../../packages/treetime/src/ancestral/fitch_sub.rs#L38-L72) <a id="cite-1"></a>[Fitch 1971](https://doi.org/10.1093/sysbio/20.4.406) [[1](#ref-1)]. That recurrence is exact for bifurcations, but it can retain non-minimum states when a node has three or more children; the arbitrary-degree recurrence retains states with minimum subtree cost <a id="cite-2"></a>[Hartigan 1973](https://doi.org/10.2307/2529676) [[2](#ref-2)].

For a root polytomy with child states $C$, $C$, and $A$, the intersection is empty and the union is $\{A,C\}$. `fn StateSet::get_one()` selects the lowest encoded state, $A$, producing two changes although choosing $C$ produces one [packages/treetime-primitives/src/bitset128.rs#L153-L159](../../packages/treetime-primitives/src/bitset128.rs#L153-L159). Sparse marginal initialization inherits the selected Fitch state.

V0 uses the same intersection-or-union recurrence [packages/legacy/treetime/treetime/treeanc.py#L638-L655](../../packages/legacy/treetime/treetime/treeanc.py#L638-L655). Correcting v1 therefore improves scientific correctness while intentionally diverging from the reference output, which requires explicit approval.

## Potential solutions

- O1. Use an exact finite-state Sankoff recurrence at nodes with three or more children and retain every minimum-cost parent state; keep the binary Fitch fast path.
- O2. Resolve multifurcations into bifurcations before Fitch reconstruction. This couples ancestral output to a topology transformation and can introduce arbitrary resolution choices.
- O3. Preserve v0 parity. This retains known non-minimum ancestral states and mutation placements.

## Recommendation

Approve O1 as an intentional v0 divergence. Validate it with the explicit $C,C,A$ counterexample, exhaustive finite-state scores on generated multifurcations, and a golden comparison that records the approved difference. No implementation ticket is ready until that scientific-output divergence is approved.

## Related issues

- [N-ancestral-fitch-site-classification-parallel-regression.md](N-ancestral-fitch-site-classification-parallel-regression.md)
- [N-test-coverage-gaps.md](N-test-coverage-gaps.md)

## Glossary

1. <a id="gloss-1"></a> **Minimum-change method.** An ancestral-state reconstruction method that minimizes the number of state transitions on a fixed tree topology. [↩](#gloss-use-1)

## References

1. <a id="ref-1"></a> Fitch, Walter M. 1971. "Toward defining the course of evolution: Minimum change for a specific tree topology." _Systematic Biology_ 20:406-416. https://doi.org/10.1093/sysbio/20.4.406 [↩](#cite-1)
2. <a id="ref-2"></a> Hartigan, John A. 1973. "Minimum mutation fits to a given tree." _Biometrics_ 29:53-65. https://doi.org/10.2307/2529676 [↩](#cite-2)
