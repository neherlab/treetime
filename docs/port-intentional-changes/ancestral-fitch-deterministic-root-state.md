# Fitch root ambiguity: deterministic selection

| Property    | Value                                                                                                                                                            |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Intentional deviation from v0                                                                                                                                    |
| v1 Location | `fitch_forward()` (`#fitch_forward`) at [`packages/treetime/src/commands/ancestral/fitch.rs#L316`](../../packages/treetime/src/commands/ancestral/fitch.rs#L316) |
| v0 Location | [`packages/legacy/treetime/treetime/treeanc.py#L606-L616`](../../packages/legacy/treetime/treetime/treeanc.py#L606-L616)                                         |
| Affects     | Root sequence at positions with tied Fitch state sets                                                                                                            |
| Datasets    | Any tree where Fitch backward pass leaves root positions ambiguous (multiple equally parsimonious states)                                                        |

## Background

Fitch parsimony (Fitch 1971) is a two-pass algorithm on phylogenetic trees.
The backward pass (leaves to root) computes the set of equally parsimonious
states at each node: intersection of children's sets if non-empty, union
otherwise. The forward pass (root to leaves) resolves ambiguity by picking the
parent's state if it is in the child's set.

At the root, there is no parent to constrain the choice. When the backward pass
leaves multiple states in the root's set (e.g. `{A, T}` are equally
parsimonious), the algorithm must pick one. This choice does not affect the
parsimony score but determines the reconstructed root sequence, which propagates
through the forward pass to all descendants.

## v0: random selection

v0 resolves root ambiguity via the per-instance RNG seeded by `rng_seed`:

```python
self.tree.root._cseq = np.array(
    [k[self.rng.integers(len(k)) if len(k) > 1 else 0] for k in self.tree.root.state]
)
```

For each ambiguous root position, `self.rng.integers(len(k))` picks a uniformly
random index from the Fitch state set. Different seeds produce different root
sequences.

## v1: deterministic first-element

v1 uses `BitSet128::get_one()` (`#get_one`), which returns the lowest-index
set bit:

```rust
for (pos, states) in variable {
    sequence[*pos] = states.get_one();
}
```

For nucleotide alphabets, this produces a systematic bias toward
alphabetically earlier states. Given tied states `{A, T}`, v1 always picks `A`
(index 0). v0 picks `A` or `T` with equal probability per the RNG seed.

## Rationale

v1 prioritizes reproducibility: identical input always produces identical output,
independent of any seed parameter. The bias toward lower-index states is
acceptable because:

- Fitch parsimony is a compression step. The marginal ML pass (which follows
  Fitch in the standard pipeline) recomputes proper posterior probabilities,
  overriding the Fitch root sequence.
- The parsimony score is unaffected by the root choice.
- For standalone Fitch usage (`--method-anc parsimony`), deterministic output
  is easier to validate and debug than seed-dependent output.

## Numerical impact

The root sequence differs at ambiguous positions. All downstream Fitch-assigned
sequences at internal nodes can differ where their state sets included the
root's chosen state. The marginal ML pass is unaffected because it computes
profiles from the GTR model independently of Fitch assignments.
