# Glossary

This glossary is optimized for TreeTime and reflects the local terminology used in this codebase and its documentation. Bioinformatics terminology is often poorly defined, with the same concept called different things by different communities and tools. This is our best effort at a clear, consistent set of definitions for the concepts that matter in TreeTime's context. We do not pretend these definitions are universally accepted or standard, and we are not trying to replace existing terminology -- just to be precise about what we mean when we use these terms.

Terms used throughout the [iterative tree refinement book](_index.md). Ordered alphabetically.

## Adaptive LASSO

An extension of [LASSO](#lasso) where each coefficient gets its own penalty weight, computed from a preliminary estimate. Branches with small initial estimates receive heavier penalties, driving them to exactly zero.

Achieves [oracle properties](#oracle-properties). Introduced by <a id="cite-22"></a>[Zou 2006](https://doi.org/10.1198/016214506000000735) [[22](#ref-22)], applied to phylogenetics by <a id="cite-21"></a>[Zhang, Dinh, and Matsen 2021](https://doi.org/10.1080/01621459.2020.1778481) [[21](#ref-21)].

## Alignment

An arrangement of biological sequences (DNA, RNA, or protein) in rows, with columns representing corresponding positions across sequences. Gaps are inserted where sequences differ in length due to insertions or deletions.

The alignment is the primary input to all tree refinement operations. Each column is an independent observation under the standard phylogenetic likelihood model.

## aLRT (Approximate Likelihood Ratio Test)

A fast branch support method that tests whether a branch has zero length by comparing the best tree against the second-best [NNI](#nni-nearest-neighbor-interchange) rearrangement.

Uses the [Self-Liang distribution](#self-liang-distribution) as the null <a id="cite-1"></a>[Anisimova and Gascuel 2006](https://doi.org/10.1080/10635150600755453) [[1](#ref-1)].

## Alternating Optimization

An iterative strategy that fixes some parameters, optimizes others, then swaps.

In phylogenetics: fix [branch lengths](#branch-length) and reconstruct ancestral states ([E-step](#e-step)), then fix ancestral states and optimize branch lengths ([M-step](#m-step)).

A form of [block coordinate descent](#block-coordinate-descent). Also called coordinate maximization or ECM.

## Among-Site Rate Variation

Not all [alignment](#alignment) positions evolve at the same rate. Conserved positions (functional constraints) evolve slowly; neutral sites evolve faster.

The discrete-gamma approximation models this by dividing the continuous gamma distribution of site rates into K categories (typically K=4) of equal probability. The site likelihood becomes a weighted average over categories. The shape parameter alpha controls rate variation: small alpha = high variation, alpha approaching infinity = uniform rates.

Can be combined with any [substitution model](#substitution-model) (e.g., GTR+G). <a id="cite-24"></a>[Yang 1994](https://doi.org/10.1007/BF00160154) [[24](#ref-24)].

## Bifurcation

An internal node in a [phylogenetic tree](#phylogenetic-tree) with exactly two children.

Most tree-building algorithms produce fully bifurcating trees.

## Block Coordinate Descent

An optimization method that partitions parameters into blocks and optimizes one block at a time while holding the others fixed.

Convergence guaranteed under regularity conditions <a id="cite-19"></a>[Tseng 2001](https://doi.org/10.1023/A:1017501703105) [[19](#ref-19)].

## Bootstrap (Phylogenetic)

A resampling method for assessing branch support.

[Alignment](#alignment) columns are resampled with replacement to create pseudoreplicate datasets, each used to build a tree. The fraction of trees containing a given clade is the bootstrap support value <a id="cite-7"></a>[Felsenstein 1985](https://doi.org/10.1111/j.1558-5646.1985.tb00420.x) [[7](#ref-7)].

## Boundary Problem

In statistical estimation, the situation where the true parameter value lies on the boundary of the allowed parameter space.

For [branch lengths](#branch-length), the boundary is zero. Standard ML estimation is biased away from boundaries, and standard [likelihood ratio test](#likelihood-ratio-test-lrt) distributions do not apply.

See [Self-Liang distribution](#self-liang-distribution).

## Branch Length

The evolutionary distance along an edge in a [phylogenetic tree](#phylogenetic-tree), measured in expected substitutions per site.

A branch length of 0.01 means approximately 1 [substitution](#substitution) per 100 alignment positions.

## Brent's Method

A derivative-free 1D optimization algorithm combining bisection, secant, and inverse quadratic interpolation.

Guaranteed to converge within a bracket <a id="cite-2"></a>[Brent 1973](https://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf) [[2](#ref-2)]. Used by v0 for [per-edge branch length optimization](#per-edge-optimization).

## Caterpillar Tree

A degenerate tree shape where each internal node has one [leaf](#leaf) child and one internal child, forming a chain. Every internal branch leads to exactly one terminal taxon.

Greedy [polytomy](#polytomy) resolution produces caterpillar-like subtrees when many children must be resolved sequentially, because each merge step peels off one pair. Stochastic [coalescent](#coalescent) resolution avoids this artifact.

## Coalescent

A population genetics model describing the genealogical history of a sample.

Lineages merge pairwise backward in time at rates determined by population size. The Kingman coalescent assumes binary mergers with exponentially distributed waiting times <a id="cite-11"></a>[Kingman 1982](https://doi.org/10.2307/3213548) [[11](#ref-11)].

Foundation for stochastic [polytomy](#polytomy) resolution.

## Compressed Branch

In v0's [polytomy](#polytomy) resolution, a child branch where `mutation_length >= clock_length` (the branch evolved at least as fast as the clock predicts).

These are less likely to be artifacts of arbitrary binary resolution.

Contrast with [stretched branch](#stretched-branch).

## Conditional Likelihood Vector

A vector of length s (number of states, 4 for nucleotides) maintained at each node u and each [alignment](#alignment) position i during the [pruning algorithm](#pruning-algorithm).

Entry `L_u(i)[a]` is the probability of observing all the data in u's subtree, given that node u has state a at position i. At [leaves](#leaf), the vector is 1 for the observed state and 0 elsewhere (or all 1s for ambiguous states).

In TreeTime, conditional likelihood vectors are called [messages](#message).

## Continuous-Time Markov Chain (CTMC)

A stochastic process on a finite set of states (e.g., {A, C, G, T}) where transitions occur continuously over time according to a [rate matrix](#rate-matrix) Q.

The Markov property: the probability of the next state depends only on the current state, not on the history. This justifies treating each branch independently when computing tree likelihoods.

All [substitution models](#substitution-model) in phylogenetics are CTMCs.

## Coordinate Sweep

One complete pass over all edges in the tree, optimizing each [branch length](#branch-length) in turn.

One sweep = one round of [per-edge optimization](#per-edge-optimization).

## Damping

Blending the new parameter value with the old one to prevent oscillation in [alternating optimization](#alternating-optimization).

The update formula:

`bl = bl_new * (1 - d^(i+1)) + bl_old * d^(i+1)`

with damping factor `d` (default 0.75) makes early iterations conservative and later iterations aggressive.

Analogous to under-relaxation in [SOR](#sor-successive-over-relaxation) <a id="cite-13"></a>[Sagulenko et al. 2018](https://doi.org/10.1093/ve/vex042) [[13](#ref-13)].

## Eigendecomposition

Factoring a matrix as `Q = V * diag(lambda) * V^-1`, where `V` contains eigenvectors and `lambda` contains eigenvalues.

For [GTR](#gtr-general-time-reversible) rate matrices, this enables efficient computation of [transition probabilities](#transition-probability-matrix):

`P(t) = V * diag(exp(lambda * t)) * V^-1`

The key property: profile-eigenvector dot products are branch-length-independent and can be precomputed once per edge.

## EM Algorithm

Expectation-Maximization: an iterative method for [maximum likelihood](#maximum-likelihood-ml) estimation when some data is missing <a id="cite-4"></a>[Dempster, Laird, and Rubin 1977](https://doi.org/10.1111/j.2517-6161.1977.tb01600.x) [[4](#ref-4)].

Alternates between computing expected values of the missing data ([E-step](#e-step)) and maximizing the likelihood given those expectations ([M-step](#m-step)). Guarantees monotonic likelihood increase.

In phylogenetics, the missing data is the ancestral sequences.

## E-Step

The expectation step of the [EM algorithm](#em-algorithm).

In phylogenetics: given current [branch lengths](#branch-length), compute the posterior probability of each ancestral state at each position at each internal node. This is [marginal reconstruction](#marginal-reconstruction).

## F81 (Felsenstein 1981)

A [substitution model](#substitution-model) with equal exchangeability rates but unequal base frequencies.

Four parameters (three free frequencies plus the overall rate). F81 sits between [JC69](#jc69-jukes-cantor) (equal frequencies) and [HKY](#hky-hasegawa-kishino-yano) (unequal frequencies plus [transition](#transitiontransversion)/[transversion](#transitiontransversion) distinction).

The per-edge log-likelihood is guaranteed [unimodal](#unimodal) under F81.

## Fitch's Algorithm

A parsimony method for ancestral state reconstruction that minimizes the total number of state changes on the tree <a id="cite-23"></a>[Fitch 1971](https://doi.org/10.2307/2412116) [[23](#ref-23)].

Two-pass dynamic programming: the backward pass (leaves to [root](#root)) computes the **state set** at each node -- the intersection of children's sets if compatible, or the union if not (incrementing the parsimony score). The forward pass (root to leaves) assigns definite states top-down.

O(n \* s) per site. With bit-set representation (v1 uses `BitSet128`), intersection and union are single AND/OR instructions.

Used in TreeTime for sequence compression, initial ancestral assignment, and mutation mapping. Not a substitute for [marginal reconstruction](#marginal-reconstruction) because it ignores [branch lengths](#branch-length) and [substitution](#substitution) rates.

## Four-Point Condition

A distance matrix can be represented as path distances in a tree if and only if for every four taxa, the two largest of the three pairwise distance sums are equal.

Characterizes when a [polytomy](#polytomy) is mathematically necessary <a id="cite-3"></a>[Buneman 1974](<https://doi.org/10.1016/0095-8956(74)90047-1>) [[3](#ref-3)].

## Grid Search

A fallback optimization method used when [Newton-Raphson](#newton-raphson) fails (second derivative non-negative).

Evaluate the likelihood at evenly spaced [branch length](#branch-length) values and pick the maximum. v1 uses 100 linearly spaced points from `0.1/L` to `1.5 * current_bl + 1/L`.

Guaranteed to find the global optimum within the grid, but slow compared to Newton-Raphson. Mitigates the risk of converging to a local optimum under multimodal likelihoods ([K2P](#k2p-kimura-2-parameter) and more complex models).

## GTR (General Time Reversible)

The most general time-reversible [substitution model](#substitution-model) for nucleotides.

Six exchangeability parameters (one per nucleotide pair) and four base frequency parameters (three free). All simpler models ([JC69](#jc69-jukes-cantor), [K2P](#k2p-kimura-2-parameter), [HKY](#hky-hasegawa-kishino-yano)) are special cases with constrained parameters <a id="cite-17"></a>[Tavaré 1986](https://doi.org/10.1090/psapm/041) [[17](#ref-17)]; <a id="cite-12"></a>[Lanave et al. 1984](https://doi.org/10.1007/BF02101990) [[12](#ref-12)].

## Hamming Distance

The number of positions at which two sequences differ.

Used as a fast initial estimate of [branch length](#branch-length): `initial_bl = hamming_distance / alignment_length`. Biased downward for long branches because multiple [substitutions](#substitution) at the same position (multiple hits) appear as a single difference.

## Hard Polytomy

A [polytomy](#polytomy) representing a genuine simultaneous divergence event: three or more descendants arose at the same time.

Rare in nature. Example: rapid radiation events.

Contrast with [soft polytomy](#soft-polytomy).

## HKY (Hasegawa-Kishino-Yano)

A [substitution model](#substitution-model) with [transition/transversion](#transitiontransversion) distinction and unequal base frequencies <a id="cite-8"></a>[Hasegawa, Kishino, and Yano 1985](https://doi.org/10.1007/BF02101694) [[8](#ref-8)].

Five parameters. The standard intermediate-complexity model.

## Information Criteria (AIC, BIC)

Model selection criteria that balance goodness of fit against model complexity.

**AIC** (Akaike Information Criterion): `AIC = -2 * log L + 2 * k`, where k is the number of parameters <a id="cite-26"></a>[Akaike 1974](https://doi.org/10.1109/TAC.1974.1100705) [[26](#ref-26)].

**BIC** (Bayesian Information Criterion): `BIC = -2 * log L + log(n) * k`, where n is the sample size <a id="cite-27"></a>[Schwarz 1978](https://doi.org/10.1214/aos/1176344136) [[27](#ref-27)].

A [polytomy](#polytomy) has one fewer parameter than the resolved tree. If the likelihood gain from resolution does not justify the extra parameter, the criterion prefers the polytomy.

## JC69 (Jukes-Cantor)

The simplest [substitution model](#substitution-model): all substitution rates equal, all base frequencies equal (1/4) <a id="cite-9"></a>[Jukes and Cantor 1969](https://doi.org/10.1016/B978-1-4832-3211-9.50009-7) [[9](#ref-9)].

One parameter. The per-edge log-likelihood is guaranteed [unimodal](#unimodal) under JC69.

## Joint Reconstruction

Finding the single most likely assignment of states to all internal nodes simultaneously, across the entire tree.

Contrast with [marginal reconstruction](#marginal-reconstruction), which optimizes each node independently. Joint reconstruction can produce different ancestral sequences at positions where the [posterior profile](#posterior-profile) is multimodal.

Implemented in v0 using the Viterbi-like algorithm of <a id="cite-25"></a>[Pupko et al. 2000](https://doi.org/10.1093/oxfordjournals.molbev.a026369) [[25](#ref-25)]. Removed in v1 as an intentional simplification.

## K2P (Kimura 2-Parameter)

A [substitution model](#substitution-model) distinguishing [transitions](#transitiontransversion) (purine-purine, pyrimidine-pyrimidine) from [transversions](#transitiontransversion) (purine-pyrimidine) <a id="cite-10"></a>[Kimura 1980](https://doi.org/10.1007/BF01731581) [[10](#ref-10)].

Two parameters. The per-edge likelihood can have multiple local maxima under K2P.

## Kappa

The [transition/transversion](#transitiontransversion) rate ratio: `kappa = alpha / beta`, where alpha is the transition rate and beta is the transversion rate.

`kappa = 1` recovers [JC69](#jc69-jukes-cantor) (no distinction). For real data, kappa ranges from ~2 (mitochondrial DNA) to ~15 (some nuclear genes). Used as the free parameter in [K2P](#k2p-kimura-2-parameter) and [HKY](#hky-hasegawa-kishino-yano).

## LASSO

Least Absolute Shrinkage and Selection Operator: an L1-regularized regression method that drives small coefficients to exactly zero, performing simultaneous estimation and variable selection <a id="cite-18"></a>[Tibshirani 1996](https://doi.org/10.1111/j.2517-6161.1996.tb02080.x) [[18](#ref-18)].

See [adaptive LASSO](#adaptive-lasso).

## Leaf

An observed sequence at the terminus of a [phylogenetic tree](#phylogenetic-tree). Also called a tip.

Each leaf represents a sampled organism: a virus isolate, a bacterial strain, a patient's tumor, or any other sequenced entity. Leaf sequences are the known data; [internal node](#root) sequences are inferred.

## Likelihood Ratio Test (LRT)

A statistical test comparing two nested models by the ratio of their maximized likelihoods.

The test statistic `2 * (log L_full - log L_restricted)` is asymptotically chi-squared distributed under the null hypothesis, except at [boundary problems](#boundary-problem) where the [Self-Liang distribution](#self-liang-distribution) applies.

## Likelihood Scaling

A numerical technique for preventing floating-point underflow in the [pruning algorithm](#pruning-algorithm) on deep trees.

Each node's [conditional likelihood vector](#conditional-likelihood-vector) is multiplied by a scaling factor to keep values in a representable range. The scaling factors are tracked and accounted for in the final log-likelihood.

Without scaling, the product of many small probabilities across hundreds of nodes underflows to zero in double-precision arithmetic.

## Long-Branch Attraction

A systematic error in parsimony-based reconstruction where convergent [substitutions](#substitution) on long branches are mistaken for shared ancestry.

Two distantly related lineages that independently acquired the same nucleotide at a position appear to share a derived state. Parsimony groups them together, producing an incorrect tree. [Maximum likelihood](#maximum-likelihood-ml) methods are less susceptible because they account for [branch lengths](#branch-length) and multiple hits.

## MAP (Maximum A Posteriori)

The single most probable value given data and model. For a probability [profile](#posterior-profile) over states at an internal node:

`MAP_u(i) = argmax_a P(state_u = a | D)`

The MAP states define the reconstructed ancestral sequence. Branch [substitutions](#substitution) are differences between the MAP states at parent and child nodes. This is how v1's `edge_subs()` determines which mutations occurred on each branch.

## Marginal Reconstruction

Computing the posterior probability distribution over ancestral states at each internal node, independently for each position.

Each node gets a probability vector (e.g., `[P(A), P(C), P(G), P(T)]`).

Contrast with [joint reconstruction](#joint-reconstruction) (finding the single most likely assignment across all nodes simultaneously).

## Maximum Likelihood (ML)

Estimating parameters by finding the values that maximize the probability of observing the data.

In phylogenetics: finding the tree [topology](#topology), [branch lengths](#branch-length), and [substitution model](#substitution-model) parameters that maximize the probability of the observed sequence [alignment](#alignment) <a id="cite-6a"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[6](#ref-6)].

## Message

TreeTime's term for the [conditional likelihood vector](#conditional-likelihood-vector) propagated during the [pruning algorithm](#pruning-algorithm).

Two messages exist at each node:

- `msg_to_parent` (backward message): computed in the leaves-to-[root](#root) pass. Carries information from the subtree below the node.
- `msg_to_child` (forward message): computed in the root-to-leaves pass. Carries information from the rest of the tree above the node.

The product of the two messages at a node gives the [posterior profile](#posterior-profile).

## M-Step

The maximization step of the [EM algorithm](#em-algorithm).

In phylogenetics: given the expected ancestral state distributions, find the [branch lengths](#branch-length) that maximize the expected complete-data log-likelihood. This is [per-edge optimization](#per-edge-optimization).

## Multifurcation

Synonym for [polytomy](#polytomy).

## Neighbor-Joining (NJ)

A distance-based tree-building algorithm that starts from a [star tree](#star-tree) (all taxa connected to a single node) and iteratively joins the closest pair, creating binary structure <a id="cite-14"></a>[Saitou and Nei 1987](https://doi.org/10.1093/oxfordjournals.molbev.a040454) [[14](#ref-14)].

Guaranteed to recover the correct tree from exact (additive) distances. O(n^3) time.

## Newton-Raphson

A root-finding/optimization method using first and second derivatives:

`x_new = x - f'(x)/f''(x)`

For [branch length](#branch-length) optimization, the analytical derivatives of the log-likelihood are available via the [eigendecomposition](#eigendecomposition).

Fast convergence near the optimum but can diverge if the second derivative is non-negative (the surface is not concave). Falls back to [grid search](#grid-search) in v1.

## NNI (Nearest Neighbor Interchange)

The simplest tree rearrangement operation: swap two subtrees adjacent to an internal edge.

Produces trees that differ by one [topology](#topology) change. Can get trapped at [polytomies](#polytomy) because neighboring topologies have near-equal likelihoods.

## one_mutation

In TreeTime, `1/L` where `L` is the [alignment](#alignment) length.

Represents the expected [branch length](#branch-length) for a single [substitution](#substitution) event. Used as a scale factor in thresholds and initial guesses.

## Oracle Properties

A statistical estimator has oracle properties if it performs as well as if the true set of zero parameters were known in advance.

Two requirements:

- Selection consistency: correctly identifies which parameters are zero.
- Asymptotic normality: estimates the non-zero parameters at the same rate as the unrestricted MLE.

## Per-Edge Optimization

Finding the optimal [branch length](#branch-length) for a single edge with fixed ancestral state distributions at both endpoints.

A 1D scalar optimization problem. The inner layer of the three-layer optimization architecture.

## Phylogenetic Tree

A branching diagram representing evolutionary relationships among organisms or sequences.

[Leaves](#leaf) represent observed sequences (samples). Internal nodes represent hypothetical ancestors. Edges (branches) connect nodes, with [branch lengths](#branch-length) representing evolutionary distance. The [root](#root) is the common ancestor of all samples.

## Polytomy

An internal node in a [phylogenetic tree](#phylogenetic-tree) with three or more children (more than two descendant branches).

A polytomy represents either genuine simultaneous divergence ([hard polytomy](#hard-polytomy)) or insufficient data to resolve the branching order ([soft polytomy](#soft-polytomy)).

Also called a [multifurcation](#multifurcation). A tree with no polytomies is fully [bifurcating](#bifurcation).

## Posterior Profile

The probability distribution over states at an internal node, given all the data in the tree.

Computed as the normalized product of the two [messages](#message) at the node:

`profile(u, i)[a] = msg_to_parent(u, i)[a] * msg_to_child(u, i)[a] / normalization`

The [MAP](#map-maximum-a-posteriori) state (argmax of the profile) defines the reconstructed ancestral state at each position.

## Pruning Algorithm

Felsenstein's dynamic programming algorithm for computing the likelihood of a tree <a id="cite-6b"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[6](#ref-6)].

Traverses from [leaves](#leaf) to [root](#root), computing [conditional likelihood vectors](#conditional-likelihood-vector) at each internal node by combining children's likelihoods through [transition probability matrices](#transition-probability-matrix).

O(n \* s^2) per site, where n is the number of nodes and s is the number of states.

## Rate Matrix (Q)

A square matrix defining the instantaneous rates of [substitution](#substitution) in a [continuous-time Markov chain](#continuous-time-markov-chain-ctmc).

Off-diagonal entry `Q_{ij}` is the instantaneous rate of change from state i to state j. Diagonal entries make each row sum to zero: `Q_{ii} = -sum_{j != i} Q_{ij}`.

The [transition probability matrix](#transition-probability-matrix) at evolutionary distance t is the matrix exponential `P(t) = exp(Q * t)`.

All [substitution models](#substitution-model) ([JC69](#jc69-jukes-cantor), [K2P](#k2p-kimura-2-parameter), [HKY](#hky-hasegawa-kishino-yano), [GTR](#gtr-general-time-reversible)) are defined by their rate matrix.

## Resolution Threshold

In [polytomy](#polytomy) resolution, the minimum likelihood gain required to merge a pair of children.

Default: 0.05 in both v0 and v1. Pairs below this threshold are left as part of the polytomy.

## Root

The single node in a [phylogenetic tree](#phylogenetic-tree) with no parent, representing the most recent common ancestor of all [leaves](#leaf).

The site likelihood at the root is computed by weighting the root's [conditional likelihood vector](#conditional-likelihood-vector) by the [stationary distribution](#stationary-distribution): `L_site(i) = sum_a pi_a * L_root(i)[a]`.

## Self-Liang Distribution

The asymptotic distribution of the [likelihood ratio test](#likelihood-ratio-test-lrt) statistic when the null hypothesis places the parameter on the boundary of the parameter space <a id="cite-15"></a>[Self and Liang 1987](https://doi.org/10.1080/01621459.1987.10478472) [[15](#ref-15)].

For a single one-sided boundary (e.g., [branch length](#branch-length) = 0), the distribution is:

`0.5 * chi^2(0) + 0.5 * chi^2(1)`

A 50-50 mixture of a point mass at zero and a chi-squared with one degree of freedom.

## Semigroup Property

For a [continuous-time Markov chain](#continuous-time-markov-chain-ctmc), the [transition probability matrices](#transition-probability-matrix) satisfy:

`P(t1 + t2) = P(t1) * P(t2)`

This matrix multiplication marginalizes over all intermediate states. It means that collapsing two consecutive edges into one requires matrix multiplication of per-edge transition matrices, not set-union of observed [substitutions](#substitution).

## Soft Polytomy

A [polytomy](#polytomy) arising from insufficient data to resolve the branching order.

The lineages did diverge at different times, but the available sequences lack the mutations needed to determine the order. Most polytomies in real data are soft. Given better data or methods, they would resolve into binary splits.

Contrast with [hard polytomy](#hard-polytomy).

## SOR (Successive Over-Relaxation)

An iterative method for solving linear systems that generalizes Gauss-Seidel by introducing a relaxation parameter omega <a id="cite-20"></a>[Young 1954](https://doi.org/10.1090/S0002-9947-1954-0059635-7) [[20](#ref-20)].

- `omega < 1`: under-relaxation (stabilizes oscillation).
- `omega > 1`: over-relaxation (accelerates convergence).

[Damping](#damping) is a form of under-relaxation.

## SPR (Subtree Pruning and Regrafting)

A tree rearrangement operation: detach a subtree and reattach it at a different point.

More powerful than [NNI](#nni-nearest-neighbor-interchange) (can reach trees that NNI cannot in one step). Used by RAxML for [topology](#topology) search.

## Star Tree

A tree where all taxa are connected to a single internal node (the [root](#root)). A complete [polytomy](#polytomy).

The star tree is the null hypothesis for polytomy testing: "no internal structure exists."

## Star Tree Paradox

[Maximum likelihood](#maximum-likelihood-ml), given enough data, always prefers a resolved ([bifurcating](#bifurcation)) tree over the true [star tree](#star-tree), because near-zero internal branches add tiny likelihood improvements at negligible parameter cost.

This means ML inherently over-resolves [polytomies](#polytomy) <a id="cite-16"></a>[Steel and Penny 2000](https://doi.org/10.1093/oxfordjournals.molbev.a026364) [[16](#ref-16)].

## Stationary Distribution

The equilibrium frequency vector `pi` of a [continuous-time Markov chain](#continuous-time-markov-chain-ctmc). As evolutionary time increases, the state frequencies converge to `pi` regardless of starting state.

For reversible [substitution models](#substitution-model), the [rate matrix](#rate-matrix) satisfies detailed balance: `pi_i * Q_{ij} = pi_j * Q_{ji}`. The [GTR](#gtr-general-time-reversible) model factors as `Q = S * diag(pi)` where S is a symmetric exchangeability matrix.

At the [root](#root) of the tree, the stationary distribution serves as the prior over root states.

## Stretched Branch

In v0's [polytomy](#polytomy) resolution, a child branch where `mutation_length < clock_length` (the branch has fewer mutations than the clock predicts).

These are candidates for being artifacts of arbitrary binary resolution by the tree builder.

Contrast with [compressed branch](#compressed-branch).

## Substitution

A change in nucleotide state at one [alignment](#alignment) position along one branch.

A to G at position 42 means the ancestral state was A and the descendant state is G.

## Substitution Model

A mathematical model describing how nucleotides change over time, formalized as a [continuous-time Markov chain](#continuous-time-markov-chain-ctmc) with [rate matrix](#rate-matrix) Q.

The model defines the probability of each nucleotide change as a function of time.

See [JC69](#jc69-jukes-cantor), [K2P](#k2p-kimura-2-parameter), [HKY](#hky-hasegawa-kishino-yano), [GTR](#gtr-general-time-reversible).

## TBR (Tree Bisection and Reconnection)

The most powerful standard tree rearrangement: bisect the tree at an edge, producing two subtrees, then reconnect them at any pair of edges.

Includes [NNI](#nni-nearest-neighbor-interchange) and [SPR](#spr-subtree-pruning-and-regrafting) as special cases.

## Topology

The branching pattern of a [phylogenetic tree](#phylogenetic-tree) -- which nodes are connected to which -- independent of [branch lengths](#branch-length).

Two trees have the same topology if they group the same sets of [leaves](#leaf) into the same clades, regardless of edge weights. Tree refinement modifies both topology (by collapsing [zero-length branches](#zero-length-branch) and resolving [polytomies](#polytomy)) and branch lengths (by [per-edge optimization](#per-edge-optimization)).

## Transition/Transversion

Two classes of nucleotide [substitution](#substitution) distinguished by the biochemistry of the change.

**Transitions** exchange purines (A and G) or pyrimidines (C and T): A-G and C-T. These require only a single tautomeric shift.

**Transversions** exchange a purine for a pyrimidine or vice versa: A-C, A-T, G-C, G-T. These require a purine-pyrimidine swap.

Transitions are 2x to 30x more frequent than transversions in real data. The ratio [kappa](#kappa) is the free parameter in [K2P](#k2p-kimura-2-parameter) and [HKY](#hky-hasegawa-kishino-yano).

Not to be confused with [transition probability matrix](#transition-probability-matrix), which refers to `P(t)`.

## Transition Probability Matrix

The matrix `P(t) = exp(Q*t)` where Q is the [rate matrix](#rate-matrix) and t is the [branch length](#branch-length).

Entry `P(t)_{ij}` is the probability that state i changes to state j over evolutionary distance t.

Satisfies the [semigroup property](#semigroup-property).

## Unimodal

Having at most one local maximum.

The per-edge log-likelihood is unimodal under [JC69](#jc69-jukes-cantor) and [F81](#f81-felsenstein-1981) models, meaning any optimization method converges to the global optimum. Under [K2P](#k2p-kimura-2-parameter) and more complex models, multiple local maxima are possible <a id="cite-5"></a>[Dinh and Matsen 2017](https://doi.org/10.1214/16-AAP1240) [[5](#ref-5)].

## Zero-Branch Penalty

In greedy [polytomy](#polytomy) resolution, a penalty for introducing a new branch with zero mutations.

Proportional to `mu * L * dt` where `mu` is the mutation rate, `L` is the [alignment](#alignment) length, and `dt` is the time span. Longer branches, higher rates, and longer alignments make zero-mutation branches less probable.

## Zero-Length Branch

A branch with [branch length](#branch-length) equal to (or indistinguishable from) zero.

Indicates either that the parent and child have identical sequences, or that the optimization converged to zero. Zero-length internal branches are candidates for collapsing (creating [polytomies](#polytomy)).

See [boundary problem](#boundary-problem), [star tree paradox](#star-tree-paradox).

## References

1. <a id="ref-1"></a> Anisimova, Maria, and Olivier Gascuel. 2006. "Approximate Likelihood-Ratio Test for Branches: A Fast, Accurate, and Powerful Alternative." _Systematic Biology_ 55(4):539-552. https://doi.org/10.1080/10635150600755453 [↩](#cite-1)
2. <a id="ref-2"></a> Brent, Richard P. 1973. _Algorithms for Minimization Without Derivatives._ Prentice-Hall. https://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf [↩](#cite-2)
3. <a id="ref-3"></a> Buneman, Peter. 1974. "A Note on the Metric Properties of Trees." _Journal of Combinatorial Theory, Series B_ 17(1):48-50. https://doi.org/10.1016/0095-8956(74)90047-1 [↩](#cite-3)
4. <a id="ref-4"></a> Dempster, Arthur P., Nan M. Laird, and Donald B. Rubin. 1977. "Maximum Likelihood from Incomplete Data via the EM Algorithm." _Journal of the Royal Statistical Society: Series B_ 39(1):1-38. https://doi.org/10.1111/j.2517-6161.1977.tb01600.x [↩](#cite-4)
5. <a id="ref-5"></a> Dinh, Vu, and Frederick A. Matsen IV. 2017. "The Shape of the One-Dimensional Phylogenetic Likelihood Function." _Annals of Applied Probability_ 27(3):1646-1677. https://doi.org/10.1214/16-AAP1240 [↩](#cite-5)
6. <a id="ref-6"></a> Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17(6):368-376. https://doi.org/10.1007/BF01734359 [↩¹](#cite-6a) [↩²](#cite-6b)
7. <a id="ref-7"></a> Felsenstein, Joseph. 1985. "Confidence Limits on Phylogenies: An Approach Using the Bootstrap." _Evolution_ 39(4):783-791. https://doi.org/10.1111/j.1558-5646.1985.tb00420.x [↩](#cite-7)
8. <a id="ref-8"></a> Hasegawa, Masami, Hirohisa Kishino, and Taka-aki Yano. 1985. "Dating of the Human-Ape Splitting by a Molecular Clock of Mitochondrial DNA." _Journal of Molecular Evolution_ 22(2):160-174. https://doi.org/10.1007/BF02101694 [↩](#cite-8)
9. <a id="ref-9"></a> Jukes, Thomas H., and Charles R. Cantor. 1969. "Evolution of Protein Molecules." In _Mammalian Protein Metabolism_, vol. 3, edited by Hans N. Munro, 21-132. Academic Press. https://doi.org/10.1016/B978-1-4832-3211-9.50009-7 [↩](#cite-9)
10. <a id="ref-10"></a> Kimura, Motoo. 1980. "A Simple Method for Estimating Evolutionary Rates of Base Substitutions Through Comparative Studies of Nucleotide Sequences." _Journal of Molecular Evolution_ 16(2):111-120. https://doi.org/10.1007/BF01731581 [↩](#cite-10)
11. <a id="ref-11"></a> Kingman, J. F. C. 1982. "The Coalescent." _Stochastic Processes and Their Applications_ 13(3):235-248. https://doi.org/10.2307/3213548 [↩](#cite-11)
12. <a id="ref-12"></a> Lanave, Cecilia, Giuliano Preparata, Cecilia Saccone, and Gabriella Serio. 1984. "A New Method for Calculating Evolutionary Substitution Rates." _Journal of Molecular Evolution_ 20(1):86-93. https://doi.org/10.1007/BF02101990 [↩](#cite-12)
13. <a id="ref-13"></a> Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-13)
14. <a id="ref-14"></a> Saitou, Naruya, and Masatoshi Nei. 1987. "The Neighbor-Joining Method: A New Method for Reconstructing Phylogenetic Trees." _Molecular Biology and Evolution_ 4(4):406-425. https://doi.org/10.1093/oxfordjournals.molbev.a040454 [↩](#cite-14)
15. <a id="ref-15"></a> Self, Steven G., and Kung-Yee Liang. 1987. "Asymptotic Properties of Maximum Likelihood Estimators and Likelihood Ratio Tests Under Nonstandard Conditions." _Journal of the American Statistical Association_ 82(398):605-610. https://doi.org/10.1080/01621459.1987.10478472 [↩](#cite-15)
16. <a id="ref-16"></a> Steel, Mike, and David Penny. 2000. "Parsimony, Likelihood, and the Role of Models in Molecular Phylogenetics." _Molecular Biology and Evolution_ 17(6):839-850. https://doi.org/10.1093/oxfordjournals.molbev.a026364 [↩](#cite-16)
17. <a id="ref-17"></a> Tavaré, Simon. 1986. "Some Probabilistic and Statistical Problems in the Analysis of DNA Sequences." _Lectures on Mathematics in the Life Sciences_ 17:57-86. https://doi.org/10.1090/psapm/041 [↩](#cite-17)
18. <a id="ref-18"></a> Tibshirani, Robert. 1996. "Regression Shrinkage and Selection via the Lasso." _Journal of the Royal Statistical Society: Series B_ 58(1):267-288. https://doi.org/10.1111/j.2517-6161.1996.tb02080.x [↩](#cite-18)
19. <a id="ref-19"></a> Tseng, Paul. 2001. "Convergence of a Block Coordinate Descent Method for Nondifferentiable Minimization." _Journal of Optimization Theory and Applications_ 109(3):475-494. https://doi.org/10.1023/A:1017501703105 [↩](#cite-19)
20. <a id="ref-20"></a> Young, David M. 1954. "Iterative Methods for Solving Partial Difference Equations of Elliptic Type." _Transactions of the American Mathematical Society_ 76(1):92-111. https://doi.org/10.1090/S0002-9947-1954-0059635-7 [↩](#cite-20)
21. <a id="ref-21"></a> Zhang, Yuxin, Vu Dinh, and Frederick A. Matsen IV. 2021. "Non-Bifurcating Phylogenetic Tree Inference via the Adaptive LASSO." _Journal of the American Statistical Association_ 116(534):858-873. https://doi.org/10.1080/01621459.2020.1778481 [↩](#cite-21)
22. <a id="ref-22"></a> Zou, Hui. 2006. "The Adaptive Lasso and Its Oracle Properties." _Journal of the American Statistical Association_ 101(476):1418-1429. https://doi.org/10.1198/016214506000000735 [↩](#cite-22)
23. <a id="ref-23"></a> Fitch, Walter M. 1971. "Toward Defining the Course of Evolution: Minimum Change for a Specific Tree Topology." _Systematic Zoology_ 20(4):406-416. https://doi.org/10.2307/2412116 [↩](#cite-23)
24. <a id="ref-24"></a> Yang, Ziheng. 1994. "Maximum Likelihood Phylogenetic Estimation from DNA Sequences with Variable Rates over Sites." _Journal of Molecular Evolution_ 39(3):306-314. https://doi.org/10.1007/BF00160154 [↩](#cite-24)
25. <a id="ref-25"></a> Pupko, Tal, Itsik Pe'er, Ron Shamir, and Dan Graur. 2000. "A Fast Algorithm for Joint Reconstruction of Ancestral Amino Acid Sequences." _Molecular Biology and Evolution_ 17(6):890-896. https://doi.org/10.1093/oxfordjournals.molbev.a026369 [↩](#cite-25)
26. <a id="ref-26"></a> Akaike, Hirotugu. 1974. "A New Look at the Statistical Model Identification." _IEEE Transactions on Automatic Control_ 19(6):716-723. https://doi.org/10.1109/TAC.1974.1100705 [↩](#cite-26)
27. <a id="ref-27"></a> Schwarz, Gideon. 1978. "Estimating the Dimension of a Model." _Annals of Statistics_ 6(2):461-464. https://doi.org/10.1214/aos/1176344136 [↩](#cite-27)
