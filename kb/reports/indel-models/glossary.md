# Glossary

[Back to index](README.md)

Consolidated glossary of domain-specific terms used across all chapters. Each term is also defined in the chapter where it first appears.

## ABC (Approximate Bayesian Computation)

<a id="gloss-abc"></a>

An inference method that bypasses likelihood computation by simulating data from prior parameter draws, comparing simulated and observed summary statistics, and accepting parameter values that produce close matches. Used by [SpartaABC](https://github.com/gilloe/SpartaABC) for SIM/RIM parameter estimation ([Loewenthal et al. 2021](https://doi.org/10.1093/molbev/msab266) [[25](references.md#ref-25)]).

First used: [Chapter 2](2-gap-treatment.md), [Chapter 5](5-empirical.md).

## Dollo principle (Dollo parsimony)

<a id="gloss-dollo"></a>

The assumption that a complex feature, once lost, cannot be regained. Applied to indels by indelMaP ([Iglhaut et al. 2024](https://doi.org/10.1093/molbev/msae109) [[26](references.md#ref-26)]): a deleted character cannot be re-inserted at the same position. This separates insertions (which create new characters) from deletions (which remove existing characters), avoiding the ambiguity of standard parsimony where a gap could be either.

First used: [Chapter 5](5-empirical.md).

## Felsenstein peeling (pruning algorithm)

<a id="gloss-peeling"></a>

The standard postorder tree traversal for computing phylogenetic likelihoods ([Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[29](references.md#ref-29)]). At each internal node, partial likelihoods from child subtrees are combined via the transition probability matrix $P(t) = e^{Qt}$. For substitution-only models, each site is processed independently. Extending peeling to indels requires tracking which alignment columns are present at each node, which is what pair HMMs and transducers accomplish ([Rivas 2008](https://doi.org/10.1371/journal.pcbi.1000172) [[9](references.md#ref-9)]).

First used: [Chapter 3](3-single-residue.md).

## Geometric distribution

<a id="gloss-geometric"></a>

$P(k) = (1-p)p^{k-1}$. The distribution implied by affine gap penalties and by TKF92 fragment lengths. Lighter-tailed than Zipf: long indels are rarer under geometric than under Zipf for the same mean. [Rivas 2005](https://doi.org/10.1186/1471-2105-6-63) [[3](references.md#ref-3)] proved that geometric instantaneous rates do not produce geometric finite-time distributions.

First used: [Chapter 2](2-gap-treatment.md).

## Gillespie algorithm

<a id="gloss-gillespie"></a>

A stochastic simulation algorithm for continuous-time Markov chains ([Gillespie 1977](https://doi.org/10.1021/j100540a008) [[30](references.md#ref-30)]). Generates exact sample paths by drawing exponentially distributed waiting times between events and choosing the next event type proportional to its rate. Used by SpartaABC to simulate indel evolution along phylogenetic branches.

First used: [Chapter 5](5-empirical.md).

## Immortal link

<a id="gloss-immortal"></a>

A TKF91 concept ([Thorne et al. 1991](https://doi.org/10.1007/BF02193625) [[7](references.md#ref-7)]): a permanent anchor at the left end of the sequence that cannot be deleted. Without it, the birth-death process could delete all characters, leaving an empty sequence with no way to restart insertions. The immortal link ensures the sequence always has at least a starting point for new insertions.

First used: [Chapter 3](3-single-residue.md).

## Pair HMM (pair hidden Markov model)

<a id="gloss-pairhmm"></a>

A hidden Markov model that generates two sequences simultaneously. States correspond to alignment column types: Match (both sequences have a residue), Insert (residue in sequence 1 only), Delete (residue in sequence 2 only). Transition probabilities encode the indel model; emission probabilities encode the substitution model. The forward algorithm computes the alignment likelihood in $O(L_1 \cdot L_2)$ time. TKF91 ([Thorne et al. 1991](https://doi.org/10.1007/BF02193625) [[7](references.md#ref-7)]), TKF92, RS07, and PIP can all be expressed as pair HMMs.

First used: [Chapter 3](3-single-residue.md).

## Phylogenetic transducer

<a id="gloss-transducer"></a>

A finite-state machine that transforms an ancestor sequence into a descendant sequence. Transducers compose: the transducer for a path through the tree is the composition of edge transducers. This enables Felsenstein-like peeling over entire sequences (not just single sites), generalizing the pruning algorithm to handle indels. Introduced by [Bradley and Holmes 2007](https://doi.org/10.1093/bioinformatics/btm402) [[28](references.md#ref-28)]; applied to indel history reconstruction by [Westesson et al. 2012](https://doi.org/10.1371/journal.pone.0034572) [[10](references.md#ref-10)].

First used: [Chapter 3](3-single-residue.md).

## POG (partial order graph)

<a id="gloss-pog"></a>

A directed acyclic graph representing a multiple sequence alignment where each path through the graph corresponds to one sequence. POGs handle insertions relative to different subsets of sequences, unlike flat matrix alignments. Used by [GRASP](https://github.com/bodenlab/GRASP) and [ProGraphMSA](https://github.com/acg-team/ProGraphMSA).

First used: [Chapter 6](6-comparison.md).

## Zipf distribution (power law)

<a id="gloss-zipf"></a>

A discrete probability distribution where the probability of observing value $k$ is $P(k) \propto k^{-a}$ for exponent $a > 1$. Larger $a$ means shorter indels dominate; smaller $a$ means long indels are more frequent. Empirical indel length data fits Zipf better than geometric in most datasets ([Cartwright 2009](https://doi.org/10.1093/molbev/msn275) [[5](references.md#ref-5)]; [Wygoda et al. 2024](https://doi.org/10.1093/bioinformatics/btae043) [[6](references.md#ref-6)]).

First used: [Chapter 2](2-gap-treatment.md), [Chapter 5](5-empirical.md).
