# Fitch recurrence is not minimum parsimony on multifurcations

The sparse Fitch backward pass implements Fitch's <a id="gloss-use-1"></a>minimum-change method <sup>[1](#gloss-1)</sup> by intersecting all child state sets and falling back to their union when the common intersection is empty [packages/treetime/src/ancestral/fitch_sub.rs#L38-L72](../../packages/treetime/src/ancestral/fitch_sub.rs#L38-L72) <a id="cite-1"></a>[Fitch 1971](https://doi.org/10.1093/sysbio/20.4.406) [[1](#ref-1)]. That recurrence is exact for bifurcations, but it can retain non-minimum states when a node has three or more children; the arbitrary-degree recurrence retains states with minimum subtree cost <a id="cite-2"></a>[Hartigan 1973](https://doi.org/10.2307/2529676) [[2](#ref-2)].

For a root polytomy with child states $C$, $C$, and $A$, the intersection is empty and the union is $\{A,C\}$. `fn StateSet::get_one()` selects the lowest encoded state, $A$, producing two changes although choosing $C$ produces one [packages/treetime-primitives/src/bitset128.rs#L153-L159](../../packages/treetime-primitives/src/bitset128.rs#L153-L159). Sparse marginal initialization inherits the selected Fitch state.

V0 uses the same intersection-or-union recurrence [packages/legacy/treetime/treetime/treeanc.py#L638-L655](../../packages/legacy/treetime/treetime/treeanc.py#L638-L655). Correcting v1 therefore improves scientific correctness while intentionally diverging from the reference output, which requires explicit approval.

## Second affected site

`resolve_variable_positions_backward` is not the only path. It iterates only positions drawn from the children's `fitch.variable` keys, which covers positions where some child carries a true IUPAC ambiguity code (`SparseNodePartition::new` populates `fitch.variable` whenever `Alphabet::is_ambiguous` holds, so leaves have entries too). Positions where children hold differing canonical fixed states and no child is variable appear in no child's `fitch.variable` and are discovered only by `resolve_fixed_positions_backward` [packages/treetime/src/ancestral/fitch_sub.rs#L86-L110](../../packages/treetime/src/ancestral/fitch_sub.rs#L86-L110), which performs no intersection at all — it accumulates a plain union of differing child states.

The two functions also both fold child states into the same `variable` map, so a position observed by both has its children counted twice. That is idempotent under a bitwise union and unsafe under any count-based rule, which makes merging them a prerequisite for the fix rather than a cleanup.

Measured on a five-leaf star with `l1..l4 = CGTT` and `l5 = AACA`: parsimony root `AACA` (16 changes) against marginal root `CGTT` (4 changes, minimum). With an internal polytomy (`(out,(l1..l5)inner)`, `out = AACA`), `inner` resolves to `AACA` at cost 4 where the minimum is 2, so the defect is not confined to the root — at non-root nodes the parent-state preference in the forward pass masks it whenever the parent happens to be right.

`get_one()` selects by lowest ASCII byte and never consults the alphabet. For `nuc` this coincides with canonical `A,C,G,T`; for `aa` it does not, since canonical order lists `*` last but `*` is ASCII 42, below `A`.

## Potential solutions

- O1. Use an exact finite-state Sankoff recurrence at nodes with three or more children and retain every minimum-cost parent state; keep the binary Fitch fast path.
- O2. Resolve multifurcations into bifurcations before Fitch reconstruction. This couples ancestral output to a topology transformation and can introduce arbitrary resolution choices.
- O3. Preserve v0 parity. This retains known non-minimum ancestral states and mutation placements.

## Recommendation

O1, approved as an intentional v0 divergence. Validate it with the explicit $C,C,A$ counterexample, exhaustive finite-state scores on generated multifurcations, and a golden comparison that records the approved difference.

O1 does not require full Sankoff cost vectors. With unit costs, `g_v(x) = C + k − count(x)` where `count(x)` is the number of children whose optimal set contains $x$ and `C = Σ_i M_i`, because a child contributes `M_i` when $x$ is in its optimal set and exactly `M_i + 1` otherwise. Hence the node's minimum-cost set is `argmax_x count(x)`, which needs only each child's optimal set plus a per-state count. The stored per-node representation remains a single `StateSet` and `SparseSeqInfo` is unchanged. The same identity shows the existing code is already exact when the intersection is non-empty (its members attain `count = k`) and when `k = 2` (an empty intersection forces `max count = 1`, so the union is the plurality set), so new behaviour is confined to `k ≥ 3` with an empty intersection.

A simple majority does not imply a unique minimum: child sets can overlap, so counts can sum above `k`. With `k` children all carrying $\{A,C\}$ both states reach `count = k`, a tie no margin over `k/2` can break. The exact test compares against `max_{y≠x} count(y)`.

Implementation is specified in [kb/tickets/ancestral-fitch-plurality-recurrence-on-multifurcations.md](../tickets/ancestral-fitch-plurality-recurrence-on-multifurcations.md).

## Related issues

- [N-ancestral-fitch-site-classification-parallel-scaling-unverified.md](N-ancestral-fitch-site-classification-parallel-scaling-unverified.md)
- [N-test-coverage-gaps.md](N-test-coverage-gaps.md)

## Glossary

1. <a id="gloss-1"></a> **Minimum-change method.** An ancestral-state reconstruction method that minimizes the number of state transitions on a fixed tree topology. [↩](#gloss-use-1)

## References

1. <a id="ref-1"></a> Fitch, Walter M. 1971. "Toward defining the course of evolution: Minimum change for a specific tree topology." _Systematic Biology_ 20:406-416. https://doi.org/10.1093/sysbio/20.4.406 [↩](#cite-1)
2. <a id="ref-2"></a> Hartigan, John A. 1973. "Minimum mutation fits to a given tree." _Biometrics_ 29:53-65. https://doi.org/10.2307/2529676 [↩](#cite-2)
