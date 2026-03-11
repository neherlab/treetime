# Optional bits in TimeTree

The simplest way to run TimeTree does not need coalescent, relaxed clock, polytomy resolution, or confidence/covariance. This simplest way is just one backward and one forward path and only needs a function that returns the likelihood of a branch having a particular length in calendar time.

The more complex iterative procedures only start surfacing once coalescent models or other complications are requested and iteration is necessary to either solve chicken/egg problems (you need the timing of nodes to calculate the coalescent likelihood) or because we do EM-style optimization where parts are optimized while others are kept fixed (e.g. optimize branches for fixed coalescent distribution, then optimize coalescent).

## Coalescent

The coalescent requires a time scale that could either be a

- constant specified by the user
- a constant that is optimized by TreeTime (opt)
- or a piecewise constant function optimized by TreeTime (skyline)

The coalescent contribution to the model is necessary for large trees with many very similar sequences.

## Relaxed clock

This feature has rarely been used in v0. It allows different branches to accumulate mutations at different rates. Technically, this is implemented by each branch having a factor that rescales the mutation rate on that branch. These factors are models by a probability distribution (either independent on each branch or coupled "auto-correlated molecular clock") and optimized as part of the iteration loop. For specific probability distributions, the optimal values can be calculated analytically.

## Polytomy resolution

This feature is important for trees with many similar sequences. In this case we often have multi-furcations. If these multi-furcations are resolved into a binary tree in a random way (as many tree builders do), the resulting random bifurcation will typically be inconsistent with the temporal ordering. TreeTime thus typically collapses zero-length branches to create polytomies and then optionally and partially resolves them in a way that is most consistent with the temporal ordering of the nodes. v0 implements two heuristics for this: One is a greedy way to always merge the pair that increases the likelihood most. This often led to a-typical tree shapes. A second method (stochastic resolve) generates more typical trees while still achieving the objective of producing clades that are compatible with the time ordering.

## Confidence and rate variation

An important feature of timetree analysis is to estimate the confidence intervals of the estimates of dates of internal nodes of the tree. At least two sources of noise contribute to this uncertainty: the randomness of mutation accumulation, and the uncertainty of the rate estimate itself. The former is well described by the gtr model and the branch length distributions and the output of the backward/forward pass contains all of this uncertainty. However, the forward/backward path is run with a fixed rate and that rate estimate itself is noisy. in v0, we deal with this by rerunning the time tree forward/backward pass with a rate that is one standard deviation higher or lower and then combining the estimates of the central and high/low rates.
Uncertainty is typically output to file as a table with upper/lower bounds for each node. In the auspice json uncertainty can be specified.
