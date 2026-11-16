# `--sample-from-profile all` draws each node independently

With `--sample-from-profile all`, the marginal reconstruction draws the state of every node from that node's own marginal posterior (`fn sample_from_profile()` [packages/treetime/src/partition/marginal/sample.rs#L39-L55](../../packages/treetime/src/partition/marginal/sample.rs#L39-L55), called with each node's profile in the dense and the sparse partition). The draws of a parent and its child are independent of each other.

A sampled set of node states is therefore not a sample of a history: the mutations between sampled nodes are not sampled mutations, a parent and child can get a pair of states that the joint posterior makes improbable, and mutation counts over many samples do not estimate the expected counts. This affects nucleotide and amino-acid runs. The `--help` text says only "sample from the posterior at every node" ([packages/app-commands/src/commands/shared/method_anc.rs#L69-L75](../../packages/app-commands/src/commands/shared/method_anc.rs#L69-L75)).

## Fix direction

Draw from the joint posterior by pre-order sampling, as BEAST does ([AncestralStateBeagleTreeLikelihood.java#L414-L530](https://github.com/beast-dev/beast-mcmc/blob/d1a0c7ef52ea1987c869e0b34889ffac62f4da8a/src/dr/evomodel/treelikelihood/AncestralStateBeagleTreeLikelihood.java#L414-L530)): the root from its posterior, then each child from $P(t)[\text{parent state}, \cdot]$ times the message from the child's subtree. With codon decoding, each draw takes one triplet per codon ([kb/proposals/codon-aware-ancestral-reconstruction.md](../proposals/codon-aware-ancestral-reconstruction.md), axis 14).

> [!IMPORTANT]
> **Investigation required.** Pre-order sampling needs the message from each child's subtree after the passes. Check which messages the dense and the sparse storage keep, and the cost of keeping them.

## Validation

On a small tree, the frequencies of sampled parent-child state pairs over many draws match the joint posterior of the pair from brute-force enumeration.
