# Dense partitions lack Fitch compression

Dense partitions skip Fitch compression entirely. `compress_sequences()` ([packages/treetime/src/commands/ancestral/fitch.rs#L526-L544](../../packages/treetime/src/commands/ancestral/fitch.rs#L526-L544)) requires the `PartitionCompressed` trait, implemented only by `PartitionFitch` and `PartitionMarginalSparse`. Dense partitions go directly to `initialize_marginal` + `update_marginal`.

Two Fitch products are absent in dense:

1. **Indels on edges.** `DenseEdgePartition.indels` ([packages/treetime/src/representation/payload/dense.rs#L37](../../packages/treetime/src/representation/payload/dense.rs#L37)) is `Vec<InDel>`, initialized empty via `Default` ([packages/treetime/src/representation/partition/marginal_dense.rs#L204](../../packages/treetime/src/representation/partition/marginal_dense.rs#L204)), never populated. The dense backward pass computes gap intersection for internal nodes ([packages/treetime/src/representation/partition/marginal_dense.rs#L221-L227](../../packages/treetime/src/representation/partition/marginal_dense.rs#L221-L227)) but does not resolve which edges gain or lose gaps.

2. **Substitutions on edges.** Fitch forward populates `fitch_subs()` on sparse edges. Dense edges have no equivalent. Dense GTR inference works around this by computing fractional mutation counts from marginal profiles, requiring a placeholder marginal pass before GTR inference (see [GTR inference chicken-and-egg problem](M-gtr-chicken-and-egg-problem.md)).

## Impact

### ~~Indel likelihood contribution is zero for dense~~ (RESOLVED)

Dense partitions now detect indels during the marginal backward/forward passes using shared `fitch_indel` logic extracted from the sparse Fitch code. `DenseEdgePartition.indels` is populated with resolved `InDel` entries, so `edge_indel_count()` returns correct values and the Poisson indel log-likelihood contributes to branch-length optimization. Tests verify dense/sparse indel count agreement.

### GTR inference requires double marginal pass

Without Fitch substitution counts, dense GTR inference ([packages/treetime/src/gtr/infer_gtr/dense.rs#L24-L38](../../packages/treetime/src/gtr/infer_gtr/dense.rs#L24-L38)) needs full marginal profiles, which need a GTR model. The workaround runs marginal with a JC69 placeholder, infers GTR from biased profiles, then reruns marginal with the real GTR. With Fitch substitution counts available, dense could infer GTR the same way sparse does: directly from integer counts, no placeholder pass needed.

## What Fitch produces and what dense needs

Fitch compression has three phases, each producing data that dense currently lacks.

### Backward pass (gap/indel detection)

`run_fitch_backward()` ([packages/treetime/src/commands/ancestral/fitch.rs#L212-L248](../../packages/treetime/src/commands/ancestral/fitch.rs#L212-L248)) resolves gap disagreements among children into `variable_indel` entries via three steps:

1. Intersect each child's gap ranges with the complement of the parent's consensus gaps to find positions where children disagree on gap presence
2. Propagate variable indel entries from children
3. If all children agree on gap, move range back to consensus gaps

This operates on `Vec<(usize, usize)>` ranges and integer counters. It does not depend on substitution tracking or sparse-specific types.

### Forward pass (indel resolution)

`run_fitch_forward()` ([packages/treetime/src/commands/ancestral/fitch.rs#L398-L439](../../packages/treetime/src/commands/ancestral/fitch.rs#L398-L439)) uses parent context to resolve variable indels into concrete `InDel::del()` / `InDel::ins()` on edges. Three cases:

1. Variable indel where deletion count + parent gap > present count: deletion edge
2. Variable indel where parent has gap but child doesn't: insertion edge
3. Consensus gaps in child not in parent (deletions) and parent gaps not in child (insertions)

### Forward pass (substitution detection)

`run_fitch_forward()` ([packages/treetime/src/commands/ancestral/fitch.rs#L340-L396](../../packages/treetime/src/commands/ancestral/fitch.rs#L340-L396)) compares Fitch-resolved child states against parent states to produce `Vec<Sub>` per edge. Dense marginal inference already reconstructs sequences from profiles (argmax), so substitutions can be derived by comparing parent and child argmax sequences after marginal reconstruction. The Fitch approach is cheaper (runs before marginal) and sufficient for GTR inference where exact posterior uncertainty is not needed.

## Design constraint: shared logic

The gap/indel arithmetic in the backward and forward passes is independent of the sparse representation. It operates on gap range lists and integer counters. Dense already computes gap intersection for internal nodes. The implementation must extract the gap-to-indel resolution logic into shared functions usable by both sparse and dense paths, not duplicate the code from `fitch.rs`.

For indels, dense does not need `InDel::seq` (the inserted sequence content) because the full sequence is available from the dense reconstruction. A `track_sequence: bool` parameter or making `InDel::seq` optional would let dense skip storing inserted sequences.

For substitutions, two approaches exist:

- Run the Fitch substitution logic over dense sequences (requires building `SparseSeqInfo`-like structures from dense data)
- Derive substitutions from argmax of dense marginal profiles after reconstruction (simpler, uses data already computed)

The substitution approach affects whether dense GTR inference can avoid the double marginal pass. Only the Fitch approach provides substitution counts before marginal inference.

## Double-counting caveat

When dense and sparse partitions represent the same alignment, `edge_indel_count()` sums across all partitions. If dense indel detection is added, partition-aware deduplication is needed to avoid doubling the Poisson indel contribution. See [optimize-indel-contribution-to-likelihood](../port-intentional-changes/optimize-indel-contribution-to-likelihood.md) double-counting caveat.

## Related

- [GTR inference chicken-and-egg problem](M-gtr-chicken-and-egg-problem.md) - the double marginal pass that Fitch subs for dense would eliminate
- [Dummy GTR initialization pattern across commands](M-core-dummy-gtr-initialization.md) - the structural manifestation of the circular dependency
- [Gap handling not implemented](H-core-gap-handling-not-implemented.md) - terminal gap filling, orthogonal but interacts with indel detection
- [Dense and sparse partition types have structural and naming asymmetries](L-representation-dense-sparse-partition-asymmetry.md) - broader asymmetry between dense and sparse
- [Indel composition missing during edge merge](M-representation-indel-composition-missing.md) - downstream indel handling
- [Indel contribution to branch length likelihood](../port-intentional-changes/optimize-indel-contribution-to-likelihood.md) - the Poisson model that reads `edge_indel_count()`
- [Dense and sparse sequence representation](../port-intentional-changes/sequence-representation-dense-sparse.md) - architectural context
