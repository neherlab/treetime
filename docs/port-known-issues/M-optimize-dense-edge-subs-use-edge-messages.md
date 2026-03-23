# Dense edge substitution counting uses edge messages instead of node posteriors

Dense `edge_subs()` counts branch mutations by comparing the argmax of `edge.msg_to_parent` and `edge.msg_to_child`. Those are partial marginal messages on the edge, not the final parent and child node posteriors. This can report mutations that are absent in the reconstructed parent-child states and can miss mutations that are present in the reconstructed states.

## Problem

The dense implementation of `edge_subs()` in [`packages/treetime/src/representation/partition/marginal_dense.rs#L67-L100`](../../packages/treetime/src/representation/partition/marginal_dense.rs#L67-L100) iterates over edge-local messages and compares their argmax states:

```rust
for (pos, parent, child) in izip!(
  0..edge.msg_to_parent.dis.nrows(),
  edge.msg_to_parent.dis.rows(),
  edge.msg_to_child.dis.rows()
) {
  let parent_state = self.alphabet.char(argmax_first(&parent).unwrap_or(0));
  let child_state = self.alphabet.char(argmax_first(&child).unwrap_or(0));
  if parent_state != child_state {
    subs.push(Sub::new(parent_state, pos, child_state)?);
  }
}
```

The same file already documents why this is suspicious: the comment at [`packages/treetime/src/representation/partition/marginal_dense.rs#L68-L72`](../../packages/treetime/src/representation/partition/marginal_dense.rs#L68-L72) states that these are partial messages, not full posteriors.

In dense marginal inference, the two quantities have different meanings:

- `edge.msg_to_parent` is the subtree-side message evaluated at the parent end of the edge. It is written during the backward pass at [`packages/treetime/src/representation/partition/marginal_dense.rs#L260-L272`](../../packages/treetime/src/representation/partition/marginal_dense.rs#L260-L272).
- `edge.msg_to_child` is the outgroup-side message for one child edge. It is written during the forward pass at [`packages/treetime/src/representation/partition/marginal_dense.rs#L303-L315`](../../packages/treetime/src/representation/partition/marginal_dense.rs#L303-L315).
- `self.nodes[&node_key].profile` is the full marginal posterior at a node, obtained after combining the incoming messages in the forward pass at [`packages/treetime/src/representation/partition/marginal_dense.rs#L277-L300`](../../packages/treetime/src/representation/partition/marginal_dense.rs#L277-L300).

Scientifically, mutation counting for a branch should compare the reconstructed state at the parent node with the reconstructed state at the child node. Comparing partial messages instead compares two incomplete views of the branch state, evaluated at different places in the message-passing graph.

## Why this is wrong

In marginal ancestral reconstruction, a node posterior is proportional to the product of independent information sources. In this implementation, the dense node posterior is formed by combining:

- the subtree-side evidence propagated upward from descendants
- the outgroup-side evidence propagated downward from the rest of the tree

The branch mutation set should therefore be derived from the final node posteriors at the two endpoints of the edge, not from one partial message from each direction.

Sparse code already follows the correct model. `PartitionMarginalSparse::edge_subs()` in [`packages/treetime/src/representation/partition/marginal_sparse.rs#L272-L295`](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L272-L295) reconstructs the current parent and child node states and compares those states directly. Its helper `node_state_at()` in [`packages/treetime/src/representation/partition/marginal_sparse.rs#L102-L165`](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L102-L165) makes the semantics explicit: branch mutations are differences between the current node states, not differences between edge-local messages.

## Example

Consider one nucleotide site on an edge:

- `edge.msg_to_parent = [0, 0, 1, 0]`, so the subtree-side message strongly supports `G`
- `edge.msg_to_child = [0.25, 0.25, 0.25, 0.25]`, so the outgroup-side message is completely uninformative

Current dense `edge_subs()` takes `argmax_first` on both vectors:

- parent-side argmax: `G`
- child-side argmax: `A`, because uniform ties resolve to the first state

The current code therefore records `G -> A` as a substitution.

That is not the reconstructed branch state. The child node posterior is formed by multiplying the child's subtree message with the outgroup message. Because the outgroup message is uniform, the child posterior remains `[0, 0, 1, 0]`, so the child state is also `G`. The correct branch state is `G -> G`, which means no substitution.

The same problem appears whenever one side of the branch is ambiguous or weakly informative. The bug is not limited to exact ties. Any case where the partial-message argmax differs from the endpoint posterior argmax can create a false branch mutation or hide a real one.

## Blast radius

### Affected production code

- `PartitionMarginalDense::edge_subs()` in [`packages/treetime/src/representation/partition/marginal_dense.rs#L67-L100`](../../packages/treetime/src/representation/partition/marginal_dense.rs#L67-L100)

### User-visible effect

On `dev`, the direct production blast radius of the exact `edge_subs()` bug is smaller than it was on `rust`. `initial_guess_mixed()` in [`packages/treetime/src/commands/optimize/optimize_unified.rs#L262-L285`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L262-L285) no longer calls `edge_subs()` for dense partitions. It now calls `edge_initial_differences()` through `PartitionOptimizeOps` in [`packages/treetime/src/commands/optimize/partition_ops.rs#L40-L44`](../../packages/treetime/src/commands/optimize/partition_ops.rs#L40-L44).

The exact defect therefore remains in production code, but no current `optimize` CLI path calls it. Exhaustive search on `dev` found no production caller that currently counts dense branch substitutions by taking argmax on `msg_to_parent` and `msg_to_child`.

The adjacent user-visible issue on `dev` is different. Dense initial branch-length seeding now uses `PartitionMarginalDense::edge_initial_differences()` in [`packages/treetime/src/representation/partition/marginal_dense.rs#L111-L135`](../../packages/treetime/src/representation/partition/marginal_dense.rs#L111-L135), which computes a soft Hamming distance from the same edge messages. That behavior is documented as intentional in [Dense initial branch length guess uses soft Hamming on per-edge messages](../port-intentional-changes/optimize-dense-initial-guess-soft-hamming.md).

That intentional change removed the old hard-argmax symptom from the `optimize` seed, but it did not answer the narrower correctness question tracked here: `edge_subs()` still reports branch mutations from partial edge messages instead of endpoint node posteriors.

Before the `dev` change, the old behavior on `rust` biased the dense contribution to the initial branch-length guess. The effect was largest on edges with ambiguous or weakly informative marginal states. The branch-length optimizer could later move away from the bad seed, but the seed still affected:

- the first optimization iterate
- convergence speed
- the chance of landing in a different local basin when the likelihood surface is flat or poorly conditioned

### Affected and related tests

- `test_dense_edge_subs_excludes_gap_positions()` in [`packages/treetime/src/commands/optimize/__tests__/test_initial_guess_gaps.rs#L214-L230`](../../packages/treetime/src/commands/optimize/__tests__/test_initial_guess_gaps.rs#L214-L230) checks gap filtering only. It does not check the message-vs-posterior distinction.
- `test_initial_guess_soft_hamming` in [`packages/treetime/src/commands/optimize/__tests__/test_initial_guess_soft_hamming.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_initial_guess_soft_hamming.rs) exercises the new `dev` initial-guess path based on `edge_initial_differences()`. These tests are relevant because they codify the edge-message interpretation that replaced the old `edge_subs()`-based seed.
- `test_sparse_edge_subs_match_reconstructed_branch_differences()` in [`packages/treetime/src/commands/ancestral/__tests__/test_marginal_sparse.rs#L499-L556`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_sparse.rs#L499-L556) already encodes the correct semantics on the sparse side.
- Existing optimize convergence tests exercise adjacent behavior, but they do not assert that dense `edge_subs()` matches reconstructed parent-child differences.
- `test_initial_guess_gtr_messages` in [`packages/treetime/src/commands/optimize/__tests__/test_initial_guess_gtr_messages.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_initial_guess_gtr_messages.rs) is also adjacent. It verifies that dense edge messages are refreshed after the dummy GTR is replaced, which matters for the new `dev` initial-guess path.

## Related code and exhaustive search

### Exact occurrence

The exact production defect exists at one site only:

- [`packages/treetime/src/representation/partition/marginal_dense.rs#L67-L100`](../../packages/treetime/src/representation/partition/marginal_dense.rs#L67-L100)

### Similar code that is correct

- Sparse branch mutation counting in [`packages/treetime/src/representation/partition/marginal_sparse.rs#L272-L295`](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L272-L295) compares reconstructed node states, not edge messages.

### Similar code that is not this bug

- Dense optimization coefficients in [`packages/treetime/src/commands/optimize/optimize_unified.rs#L52-L65`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L52-L65) and [`packages/treetime/src/commands/optimize/optimize_dense.rs#L67-L75`](../../packages/treetime/src/commands/optimize/optimize_dense.rs#L67-L75) also consume `msg_to_parent` and `msg_to_child`, but that is scientifically correct there. Branch likelihood is defined on the edge messages. The bug is specific to using those same messages as if they were endpoint states for mutation counting.
- Dense initial branch-length seeding in [`packages/treetime/src/representation/partition/marginal_dense.rs#L111-L135`](../../packages/treetime/src/representation/partition/marginal_dense.rs#L111-L135) is a close sibling concern on `dev`. It also consumes `msg_to_parent` and `msg_to_child`, but it computes a soft overlap score rather than hard branch substitutions. This is documented as an intentional change, not this exact bug. If the team later decides that initial branch differences must also be evaluated at a common node posterior, that intentional-change entry should be revisited.
- Dense GTR inference in [`packages/treetime/src/gtr/infer_gtr/dense.rs#L49-L80`](../../packages/treetime/src/gtr/infer_gtr/dense.rs#L49-L80) also consumes edge messages, but it computes a normalized joint distribution over parent and child states and sums expected counts. That code integrates over uncertainty. It does not convert partial messages into discrete branch substitutions, so it is not this bug.

### Related known issues

- [Branch mutations have no unified API across partition types](M-core-branch-mutations-no-unified-api.md)

The current issue is narrower and more concrete than the adjacent intentional-change entry. Even if v0 did not exist, comparing partial edge messages instead of node posteriors would still be wrong for mutation counting.

## Proposed solutions

### S1. Compare parent and child node posteriors directly

Change dense `edge_subs()` to read:

- `self.nodes[&parent_key].profile.dis`
- `self.nodes[&child_key].profile.dis`

Then compare the argmax state at each non-gap position.

Why this works:

- it matches the scientific meaning of a branch mutation
- it matches the node posteriors used for dense sequence reconstruction
- it aligns dense semantics with the existing sparse implementation

Pros:

- minimal conceptual change
- no extra reconstruction pass
- direct reuse of already-computed node posteriors
- easy to verify against reconstructed sequences

Cons:

- only valid after the forward pass has populated node posteriors
- still uses hard MAP state differences, not fractional expected mutation counts

### S2. Reconstruct parent and child sequences, then diff them

Reuse the same sequence assignment logic used by dense ancestral reconstruction, then compare the parent and child sequences position by position.

Why this works:

- it guarantees consistency with the sequences that users and downstream code actually see

Pros:

- strongest semantic consistency with reconstructed sequence output
- easy to reason about in tests

Cons:

- more allocation and work per edge
- duplicates sequence reconstruction effort when only a branch diff is needed
- less direct than reading node posteriors

### Recommended direction

S1 is the better fix. The production bug is that `edge_subs()` reads the wrong data structure. Switching to endpoint node posteriors corrects that mistake directly and keeps the implementation small.

## Implementation notes

- Do not change dense likelihood coefficient code in [`packages/treetime/src/commands/optimize/optimize_unified.rs#L52-L65`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L52-L65) or [`packages/treetime/src/commands/optimize/optimize_dense.rs#L67-L75`](../../packages/treetime/src/commands/optimize/optimize_dense.rs#L67-L75) as part of this fix. Those paths use edge messages for likelihood integration, which is correct.
- Do not change dense GTR inference in [`packages/treetime/src/gtr/infer_gtr/dense.rs#L49-L80`](../../packages/treetime/src/gtr/infer_gtr/dense.rs#L49-L80) as part of this fix. That path already computes expected joint counts rather than MAP state differences.
- Keep `edge_subs()` and `edge_initial_differences()` as separate operations unless the team explicitly decides to revisit the `dev` intentional change. They answer different scientific questions and should not silently collapse into one abstraction.

## Proposed tests and verification

### T1. Dense regression test for ambiguous outgroup message

Construct a small dense partition where:

- the parent node posterior is sharp at `G`
- the child node posterior is also sharp at `G`
- the child edge's `msg_to_child` is uniform at the same site

Expected result:

- current code would report one substitution
- fixed code must report zero substitutions

### T2. Dense parity test against reconstructed parent-child sequences

For each edge after `update_marginal()`:

- derive `edge_subs()`
- reconstruct the parent and child MAP sequences from node posteriors
- compare `edge_subs()` against the explicit parent-child sequence diff

This should mirror the sparse test pattern in [`packages/treetime/src/commands/ancestral/__tests__/test_marginal_sparse.rs#L499-L556`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_sparse.rs#L499-L556).

### T3. Initial-guess regression

For the exact `edge_subs()` bug, create an edge where the dense child-side partial message is ambiguous but the final child node posterior agrees with the parent. Assert that `edge_subs()` reports zero substitutions.

If the team decides to revisit the `dev` soft-Hamming design, add a separate test that compares `edge_initial_differences()` against a node-posterior-based overlap measure and documents the intended relationship between the two.

### T4. Existing gap tests remain green

Keep `test_dense_edge_subs_excludes_gap_positions()` in [`packages/treetime/src/commands/optimize/__tests__/test_initial_guess_gaps.rs#L214-L230`](../../packages/treetime/src/commands/optimize/__tests__/test_initial_guess_gaps.rs#L214-L230) to verify that the fix does not regress gap filtering.
