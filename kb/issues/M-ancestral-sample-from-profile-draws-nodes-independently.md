# `--sample-from-profile all` draws each node independently

With `--sample-from-profile all`, the marginal reconstruction draws the state of every node from that node's own marginal posterior (`fn sample_from_profile()` [packages/treetime/src/partition/marginal/sample.rs#L39-L55](../../packages/treetime/src/partition/marginal/sample.rs#L39-L55), called through `fn resolve_profile()` with each node's profile by `fn sample_sequences()` of the dense partition [packages/treetime/src/partition/marginal/dense/partition.rs#L181-L191](../../packages/treetime/src/partition/marginal/dense/partition.rs#L181-L191) and of the sparse partition [packages/treetime/src/partition/marginal/sparse/partition.rs#L246](../../packages/treetime/src/partition/marginal/sparse/partition.rs#L246)). The draws of a parent and its child are independent of each other.

A sampled set of node states is therefore not a sample of a history: the mutations between sampled nodes are not sampled mutations, a parent and child can get a pair of states that the joint posterior makes improbable, and mutation counts over many samples do not estimate the expected counts. This affects nucleotide and amino-acid runs. The `--help` text says only "sample from the posterior at every node" ([packages/app-commands/src/commands/shared/method_anc.rs#L93](../../packages/app-commands/src/commands/shared/method_anc.rs#L93)).

## Fix direction

Draw from the joint posterior by pre-order sampling, as BEAST does ([AncestralStateBeagleTreeLikelihood.java#L414-L530](https://github.com/beast-dev/beast-mcmc/blob/d1a0c7ef52ea1987c869e0b34889ffac62f4da8a/src/dr/evomodel/treelikelihood/AncestralStateBeagleTreeLikelihood.java#L414-L530)): the root from its posterior, then each child from $P(t)[\text{parent state}, \cdot]$ times the message from the child's subtree. With codon decoding, each draw takes one triplet per codon ([kb/proposals/codon-aware-ancestral-reconstruction.md](../proposals/codon-aware-ancestral-reconstruction.md), axis 14).

> [!IMPORTANT]
> **Investigation required.** Pre-order sampling needs the message from each child's subtree after the passes. Check which messages the dense and the sparse storage keep, and the cost of keeping them.

## Validation

On a small tree, the frequencies of sampled parent-child state pairs over many draws match the joint posterior of the pair from brute-force enumeration.

## Proposed implementation plan

1. Resolve the **Investigation required** block first: list which per-node messages the dense (`packages/treetime/src/partition/marginal/dense/`) and the sparse (`packages/treetime/src/partition/marginal/sparse/`) marginal passes keep after the backward pass (the subtree message of each child, per site), and measure the extra memory of keeping them on a large dataset (`data/sc2/4500`). The steps below assume the messages are kept or recomputed per child.
2. Add a pre-order sampler beside `fn sample_from_profile()` in `packages/treetime/src/partition/marginal/sample.rs`: draw the root from its posterior; for each child in pre-order, draw from the normalized product of row `P(t)[parent state, :]` and the child's subtree message. Leaves with observed states keep their data.
3. Route `SampleMode::All` through the new sampler in `fn sample_sequences()` of the dense and the sparse partition, so that the internal states come from one joint draw. `SampleMode::Root` keeps its behavior (one draw at the root, argmax elsewhere). In the sparse partition, sample the variable sites with the sampler and fill the fixed and globally missing sites (`sparse/globally_missing.rs`) with the same conditional rule.
4. Keep the codon case of [kb/proposals/codon-aware-ancestral-reconstruction.md](../proposals/codon-aware-ancestral-reconstruction.md) (axis 14) out of scope until that proposal is decided.
5. Update the `--help` doc comment at `packages/app-commands/src/commands/shared/method_anc.rs#L93` to say that 'all' draws one history from the joint posterior, then regenerate the generated files (`just gen`).

### Tests

- Unit (statistical), new file `packages/treetime/src/partition/marginal/__tests__/test_sample_joint.rs`: a three-leaf tree with fixed branch lengths and JC69, seeded `StdRng`, N = 100000 draws; compare the frequencies of sampled parent-child state pairs with the joint posterior of the pair from brute-force enumeration over all internal states (oracle: the enumeration, independent of the passes). Use a chi-square bound with a stated false-fail rate.
- The same test on the dense and on the sparse partition, as an `rstest` case per representation.
- Property: for a tree where all leaves are identical and branches short, every sampled history has zero mutations with probability near one (per-run invariant: sampled states are always states with nonzero joint probability).

### Validation

`just check`, `just test-rs`, `just generated-check`; `just r treetime ancestral --sample-from-profile all` on `data/flu/h3n2/20` to check that the run completes and the mutation count per sample does not exceed the count of independent per-node draws.
