# Ancestral Reconstruction Algorithms

[Back to index](../index.md)

## Fitch Parsimony

| Property    | Value                                                                                                                                                                      |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                                                                 |
| v1 Location | [`packages/treetime/src/commands/ancestral/fitch.rs#L85-L513`](../../../packages/treetime/src/commands/ancestral/fitch.rs#L85-L513)                                        |
| v0 Location | [`packages/legacy/treetime/treetime/treeanc.py#L575-L686`](../../../packages/legacy/treetime/treetime/treeanc.py#L575-L686)                                                |
| Functions   | `fitch_backward()` (`#fitch_backward`), `fitch_forward()` (`#fitch_forward`), `run_fitch_backward()` (`#run_fitch_backward`), `run_fitch_forward()` (`#run_fitch_forward`) |
| Reference   | Fitch, W.M. (1971). "Toward Defining the Course of Evolution." Systematic Zoology, 20(4):406-416                                                                           |
| Paper URL   | https://doi.org/10.2307/2412116                                                                                                                                            |

Two-pass DP: backward (leaves-to-root) computes intersection/union of child state sets; forward (root-to-leaves) picks parent state if in child set else picks arbitrarily.

**v1 Extensions**: Sparse representation (only variable positions stored), indel handling with majority rule, `BitSet128` state sets for O(1) set operations, parallel BFS traversal.

**Modifications**: v1 uses deterministic `get_one()` (`#get_one`) for root state selection; v0 uses random selection.

---

## Marginal ML

| Property    | Value                                                                                                                                                                                      |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Type        | Well-known                                                                                                                                                                                 |
| v1 Dense    | [`packages/treetime/src/representation/partition/marginal_dense.rs#L87-L281`](../../../packages/treetime/src/representation/partition/marginal_dense.rs#L87-L281)                          |
| v1 Sparse   | [`packages/treetime/src/representation/partition/marginal_passes.rs#L16-L250`](../../../packages/treetime/src/representation/partition/marginal_passes.rs#L16-L250)                        |
| v0 Location | [`packages/legacy/treetime/treetime/treeanc.py#L762-L927`](../../../packages/legacy/treetime/treetime/treeanc.py#L762-L927)                                                                |
| Functions   | `process_node_backward()` (`#process_node_backward`), `process_node_forward()` (`#process_node_forward`), `combine_messages()` (`#combine_messages`), `propagate_raw()` (`#propagate_raw`) |
| Reference   | Felsenstein, J. (1981). "Evolutionary trees from DNA sequences." J Mol Evol, 17(6):368-376                                                                                                 |
| Paper URL   | https://doi.org/10.1007/BF01734359                                                                                                                                                         |

Belief propagation / sum-product algorithm on trees. Backward pass computes partial likelihoods via GTR matrix multiplication; forward pass computes outgroup messages via cavity/division.

**v1 vs v0**: v1 uses plain probability space; v0 uses neg-log space. v1 dense uses deterministic `argmax_first()` (`#argmax_first`) for sequence extraction; v0 uses `np.argmax()` which has undefined tie-breaking.

---

## Joint ML (Unimplemented)

See [unimplemented](../unimplemented/index.md#joint-ml) for full details.

| Property    | Value                                                                                                                                         |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                                    |
| v1 Status   | `unimplemented!()` at [`packages/treetime/src/commands/ancestral/run.rs#L194`](../../../packages/treetime/src/commands/ancestral/run.rs#L194) |
| v0 Location | [`packages/legacy/treetime/treetime/treeanc.py#L934-L1080`](../../../packages/legacy/treetime/treetime/treeanc.py#L934-L1080)                 |
| Reference   | Pupko et al. (2000). "A fast algorithm for joint reconstruction." Mol Biol Evol, 17(6):890-896                                                |

Uses traceback pointers (argmax) instead of marginalization.

---

## File Index

| File                                                                                                                                                | Algorithms                                                                       |
| --------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------- |
| [`packages/treetime/src/commands/ancestral/fitch.rs`](../../../packages/treetime/src/commands/ancestral/fitch.rs)                                   | Fitch parsimony (backward, forward, cleanup)                                     |
| [`packages/treetime/src/commands/ancestral/marginal.rs`](../../../packages/treetime/src/commands/ancestral/marginal.rs)                             | Marginal ML orchestration                                                        |
| [`packages/treetime/src/commands/ancestral/run.rs`](../../../packages/treetime/src/commands/ancestral/run.rs)                                       | Ancestral command entry point, method dispatch                                   |
| [`packages/treetime/src/representation/partition/marginal_dense.rs`](../../../packages/treetime/src/representation/partition/marginal_dense.rs)     | Dense marginal (Felsenstein pruning)                                             |
| [`packages/treetime/src/representation/partition/marginal_sparse.rs`](../../../packages/treetime/src/representation/partition/marginal_sparse.rs)   | Sparse marginal                                                                  |
| [`packages/treetime/src/representation/partition/marginal_passes.rs`](../../../packages/treetime/src/representation/partition/marginal_passes.rs)   | Sparse message passing                                                           |
| [`packages/treetime/src/representation/partition/marginal_helpers.rs`](../../../packages/treetime/src/representation/partition/marginal_helpers.rs) | `combine_messages()` (`#combine_messages`), `propagate_raw()` (`#propagate_raw`) |
