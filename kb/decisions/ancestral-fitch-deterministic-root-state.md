# Fitch Root Ambiguity: Deterministic Selection

This document describes an intentional deviation from v0 behavior in the v1 Rust implementation. The change affects root sequence reconstruction at positions where Fitch parsimony produces multiple equally optimal states. v1 resolves this in `fitch_forward()` at `packages/treetime/src/commands/ancestral/fitch.rs:315-317:`, while v0 handles it in `_fitch_reconstruction()` at `packages/legacy/treetime/treetime/treeanc.py:615-617:`. The change applies to any phylogenetic tree where the backward pass leaves root positions ambiguous.

## Fitch Parsimony Algorithm

Walter M. Fitch introduced the maximum parsimony algorithm in 1971 (Fitch, W. M. "Toward defining the course of evolution: minimum change for a specific tree topology." _Systematic Zoology_ 20(4):406-416, 1971). The algorithm finds ancestral character states that minimize the total number of state changes across a phylogenetic tree.

### The Two-Pass Algorithm

Fitch parsimony operates in two passes over the tree:

**Backward pass (leaves to root)**: At each internal node, compute the set of states that could be ancestral without requiring an extra mutation. For a node with children having state sets S1 and S2:

- If S1 and S2 overlap, the node's state set is their intersection (S1 ∩ S2)
- If S1 and S2 are disjoint, the node's state set is their union (S1 ∪ S2), and the parsimony score increases by 1

**Forward pass (root to leaves)**: Assign concrete states to each node. For non-root nodes:

- If the parent's assigned state is in this node's state set, use the parent's state (no mutation)
- Otherwise, pick any state from this node's state set (one mutation)

### Root Ambiguity

At the root, there is no parent to constrain the choice. When the backward pass produces a state set with multiple elements (e.g., {A, T}), any choice is equally parsimonious. Mathematically, if the root state set R has |R| > 1 elements, all states in R produce the same minimum parsimony score. The algorithm must select one state to proceed with the forward pass.

This selection does not affect the parsimony score, but it determines the reconstructed root sequence. Since the forward pass propagates the root's state downward (preferring parent states when possible), the root choice can affect ancestral sequences at all internal nodes.

## v0 Implementation: Random Selection

v0 resolves root ambiguity using a random number generator in `packages/legacy/treetime/treetime/treeanc.py:615-617:`:

```python
self.tree.root._cseq = np.array(
    [k[self.rng.integers(len(k)) if len(k) > 1 else 0] for k in self.tree.root.state]
)
```

For each position, `self.tree.root.state[i]` contains a list of equally parsimonious states. When this list has more than one element, `self.rng.integers(len(k))` selects a uniformly random index. The RNG is initialized at construction time in `packages/legacy/treetime/treetime/treeanc.py:163:` with `np.random.default_rng(seed=rng_seed)`.

Different seeds produce different root sequences, and thus different ancestral reconstructions at internal nodes.

Note: The v0 logging code at `packages/legacy/treetime/treetime/treeanc.py:612:` claims to choose `state[amb][0]` but the actual selection at line 616 uses the RNG. This is a logging inconsistency in v0.

## v1 Implementation: Deterministic Selection

v1 resolves root ambiguity deterministically in `packages/treetime/src/commands/ancestral/fitch.rs:315-317:`:

```rust
for (pos, states) in variable {
    sequence[*pos] = states.get_one();
}
```

The `states` variable is a `BitSet128`, which stores character states as bits in a 128-bit integer where bit position equals the ASCII code of the character. The `get_one()` method at `packages/treetime-primitives/src/bitset128.rs:165-167:` returns the first set bit:

```rust
pub fn get_one(&self) -> AsciiChar {
    self.get_one_maybe().expect("BitSet128 is empty")
}
```

This calls `first()` at `packages/treetime-primitives/src/bitset128.rs:153-155:`:

```rust
pub fn first(&self) -> Option<AsciiChar> {
    (!self.is_empty()).then_some(AsciiChar::from_byte_unchecked(self.bits.trailing_zeros() as u8))
}
```

The `trailing_zeros()` operation returns the index of the lowest set bit, which corresponds to the character with the smallest ASCII code. For DNA nucleotides, the ASCII codes are: A=65, C=67, G=71, T=84. Given an ambiguous state set {A, T}, v1 always selects A (lowest ASCII code). Given {C, G}, v1 always selects C.

## Rationale

v1 prioritizes reproducibility: identical input always produces identical output, independent of any seed parameter. The systematic preference for lower ASCII codes is acceptable because:

1. **Fitch is a preprocessing step**: In the standard TreeTime pipeline, Fitch parsimony provides initial ancestral estimates that the marginal maximum likelihood pass subsequently refines. The ML pass computes posterior probabilities from the GTR substitution model independently of Fitch assignments.

2. **Parsimony score is unaffected**: All states in an ambiguous set are equally parsimonious by definition. The choice affects only which equally-optimal reconstruction is returned.

3. **Determinism aids validation**: For standalone parsimony usage (`--method-anc parsimony`), deterministic output simplifies testing and debugging. Comparing outputs across runs or implementations does not require controlling for random seed effects.

## Numerical Impact

The root sequence differs from v0 at ambiguous positions. All downstream Fitch-assigned sequences at internal nodes can differ where their state sets included the root's chosen state. The marginal ML pass remains unaffected because it computes profiles from the GTR model independently of Fitch assignments.
