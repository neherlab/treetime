# Joint ML ancestral reconstruction removed

v1 removes support for joint maximum likelihood ancestral reconstruction. The enum variant `MethodAncestral::Joint` (`#MethodAncestral`) in `packages/treetime/src/commands/ancestral/args.rs:8-16:` is preserved as `#[default]` for CLI compatibility with v0, but selecting it triggers `unimplemented!()` at `packages/treetime/src/commands/ancestral/run.rs:201-203:`. The v0 implementation lives in `_ml_anc_joint()` (`#_ml_anc_joint`) at `packages/legacy/treetime/treetime/treeanc.py:934-1080:`.

Currently only the `ancestral` command panics on joint. The `timetree` and `clock` commands have `method_anc` fields in their argument structs but do not dispatch on them yet.

## Background: marginal vs joint reconstruction

Ancestral sequence reconstruction infers the sequences at internal nodes of a phylogenetic tree given observed sequences at the leaves. Two maximum likelihood approaches exist: marginal and joint.

**Marginal reconstruction** computes the posterior probability distribution over states at each position independently for each internal node. At position i of node n, it asks: "Given the observed data and the substitution model, what is the probability that this position is A, C, G, or T?" The algorithm integrates (marginalizes) over all possible states at all other nodes. The result is a probability profile at each node - a vector of probabilities summing to 1 for each alignment position.

**Joint reconstruction** finds the single most probable assignment of states across all nodes simultaneously. It asks: "What configuration of ancestral sequences across the entire tree maximizes the probability of observing the leaf sequences?" The algorithm uses dynamic programming with backpointers (similar to the Viterbi algorithm for hidden Markov models) to find this global optimum. The result is a single sequence at each node - a point estimate with no uncertainty information.

The distinction matters because the most probable state at each position (from marginal) may differ from the state assigned to that position in the globally most probable configuration (from joint). When uncertainty is high, marginal reconstruction reveals this through diffuse probability profiles; joint reconstruction hides it by committing to point estimates.

## The Pupko et al. 2000 algorithm

Joint reconstruction was formalized by <a id="cite-1"></a>[Pupko et al. 2000](https://doi.org/10.1093/oxfordjournals.molbev.a026369) [[1](#ref-1)]. Their algorithm runs in O(n _ L _ k^2) time where n is the number of taxa, L is sequence length, and k is alphabet size.

The algorithm has two passes over the tree:

**Postorder pass** (leaves to root): For each node, compute two arrays indexed by alignment position and possible parent state:

- `Lx[pos, parent_state]`: Log-likelihood of the best subtree assignment when the parent has this state
- `Cx[pos, parent_state]`: Which child state achieves this maximum (backpointer)

At each position, for each possible parent state j, the algorithm finds the child state i that maximizes `log P(j -> i) + Lx_children[pos, i]`, where P(j -> i) is the transition probability from the GTR substitution model.

**Root selection**: The root state at each position maximizes `Lx[pos, state] + log(pi[state])`, where pi is the equilibrium frequency.

**Preorder pass** (root to leaves): Starting from the chosen root states, trace back through the `Cx` backpointers to recover the optimal state at each internal node. At each node, the state is `Cx[pos, parent_state]` where `parent_state` is already determined from the node's parent.

The v0 implementation at `packages/legacy/treetime/treetime/treeanc.py:960-1063:` follows this structure exactly, operating entirely in log-space to avoid numerical underflow.

## Why joint reconstruction is statistically inconsistent

<a id="cite-2"></a>[Mossel, Roch, and Steel 2009](https://doi.org/10.1109/TCBB.2008.107) [[2](#ref-2)] proved that joint ancestral maximum likelihood is statistically inconsistent. This means that even with infinite sequence data, the method can produce incorrect results.

The problem is branch length shrinkage. When jointly optimizing over ancestral states and branch lengths, the optimization systematically underestimates branch lengths. The estimator prefers shorter branches because they make fewer state changes more likely, and when ancestral states are also being optimized, the algorithm can "explain away" the evidence for longer branches by choosing convenient ancestral states.

For certain tree topologies - particularly when one pair of sister branches is long and internal branches are short - the joint method estimates internal branch lengths as exactly zero even with infinite data. This collapses internal resolution, turning a resolved tree into a star tree.

<a id="cite-3"></a>[Shaw, Dinh, and Matsen 2019](https://doi.org/10.1093/molbev/msz128) [[3](#ref-3)] strengthened this result: the only parameter values for which joint inference produces correct branch lengths lie in a set of measure zero. The bias is systematic and downward.

Marginal reconstruction avoids this problem by integrating over ancestral states rather than optimizing them. The classical Felsenstein pruning algorithm computes marginal likelihoods that are statistically consistent for branch length and topology inference.

## Why v1 uses marginal only

v1's architecture is built around probability profiles that flow through the tree during inference. The partition system stores probability distributions at each node:

- `PartitionMarginalDense`: Full probability vectors (k values per position)
- `PartitionMarginalSparse`: Probability vectors only at variable positions

Downstream algorithms depend on these probability profiles:

- GTR model inference uses expected state frequencies from profiles
- Timetree optimization monitors convergence via profile changes
- Uncertainty quantification reads probability distributions directly

Joint reconstruction produces sequences (single states per position), not profiles. Integrating it would require either:

1. Converting point estimates back to degenerate profiles (losing all uncertainty information)
2. Maintaining a parallel code path that bypasses profile-based algorithms

Neither option aligns with v1's design. The marginal approach provides strictly more information (full posteriors vs point estimates) and is statistically consistent. Joint reconstruction offers no compensating advantage for TreeTime's use cases.

## Practical impact

Running `treetime ancestral` without `--method-anc` panics because `Joint` is the default. Users must specify `--method-anc marginal` or `--method-anc parsimony`.

The `timetree` and `clock` commands accept `--method-anc` but ignore it - they always use marginal reconstruction internally. This may change in future versions.

No golden master tests exist for joint reconstruction in v1. Marginal reconstruction is the recommended method for all use cases.

## References

- <a id="ref-1"></a>Pupko, Tal, Itsik Pe'er, Ron Shamir, and Dan Graur. 2000. "A Fast Algorithm for Joint Reconstruction of Ancestral Amino Acid Sequences." _Molecular Biology and Evolution_ 17(6):890-896. https://doi.org/10.1093/oxfordjournals.molbev.a026369 [↩](#cite-1)
- <a id="ref-2"></a>Mossel, Elchanan, Sebastien Roch, and Mike Steel. 2009. "Shrinkage Effect in Ancestral Maximum Likelihood." _IEEE/ACM Transactions on Computational Biology and Bioinformatics_ 6(1):126-133. https://doi.org/10.1109/TCBB.2008.107 [↩](#cite-2)
- <a id="ref-3"></a>Shaw, David, Vu Dinh, and Frederick A. Matsen IV. 2019. "Joint Maximum Likelihood of Phylogeny and Ancestral States Is Not Consistent." _Molecular Biology and Evolution_ 36(11):2613-2619. https://doi.org/10.1093/molbev/msz128 [↩](#cite-3)
