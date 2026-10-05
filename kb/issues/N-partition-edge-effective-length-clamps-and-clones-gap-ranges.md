# `edge_effective_length()` clamps to zero and clones gap ranges

`edge_effective_length()` returns the number of alignment positions on an edge where both the parent and the child have a character. The dense and sparse implementations have two defects: they clamp an impossible negative length to zero without reporting it, and they clone both non-character range vectors on every call.

## Instances

### saturating_sub clamps to 0

- Sparse `fn edge_effective_length()` [packages/treetime/src/partition/marginal/sparse/partition.rs#L71-L82](../../packages/treetime/src/partition/marginal/sparse/partition.rs#L71-L82), clamp at L81
- Dense `fn edge_effective_length()` [packages/treetime/src/partition/marginal/dense/partition.rs#L126-L142](../../packages/treetime/src/partition/marginal/dense/partition.rs#L126-L142), clamp at L141

`self.length.saturating_sub(non_char_positions)` returns 0 when the union of non-character ranges covers more positions than the sequence length. That state breaks the invariant that the ranges lie inside the sequence, but the caller cannot tell it apart from an edge whose positions are all gaps or unknown.

The only consumer, `fn initial_guess_mixed()` [packages/treetime/src/optimize/dispatch.rs#L203-L220](../../packages/treetime/src/optimize/dispatch.rs#L203-L220), guards the division `sub_count / effective_length` with `effective_length > 0`, so no division by zero occurs. For a zero length it falls back to `one_mutation` (edges with indels) or `0.1 * one_mutation`, so a broken range invariant silently produces a plausible initial branch length. `fn gather_edge_effective_lengths()` [packages/treetime/src/optimize/gather.rs#L86](../../packages/treetime/src/optimize/gather.rs#L86) collects the values for every edge.

### range_union clones gap vectors

- Sparse call [packages/treetime/src/partition/marginal/sparse/partition.rs#L76](../../packages/treetime/src/partition/marginal/sparse/partition.rs#L76)
- Dense call [packages/treetime/src/partition/marginal/dense/partition.rs#L136](../../packages/treetime/src/partition/marginal/dense/partition.rs#L136)

`range_union(&[parent_non_char.clone(), child_non_char.clone()])` clones both range vectors because `fn range_union()` takes a slice of owned vectors. `fn range_union_iter()` [packages/treetime-utils/src/interval/range_union.rs#L9-L15](../../packages/treetime-utils/src/interval/range_union.rs#L9-L15) takes an iterator of references and needs no clone. The function runs once per edge before optimization, so the allocation pressure grows on heavily gapped alignments.

## Proposed solution

- Return an error when the non-character positions exceed the sequence length, instead of clamping
- Compute the union through `range_union_iter()` over references to the two range vectors
