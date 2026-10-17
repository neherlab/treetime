# Optional components in timetree inference

The minimal timetree pipeline runs a single backward and forward pass, requiring only a function that returns the likelihood of a branch having a particular length in calendar time. This path does not need coalescent models, relaxed clock, polytomy resolution, or confidence/covariance.

The more complex iterative procedures surface when coalescent models or other options are requested. Iteration is necessary to solve chicken-and-egg problems (node timing is needed to calculate the coalescent likelihood, but coalescent affects timing) or because the algorithm uses EM-style optimization where parts are optimized while others are kept fixed (e.g., optimize branches for fixed coalescent distribution, then optimize coalescent).

## Coalescent

The coalescent requires a time scale that is either:

- a constant specified by the user
- a constant optimized by TreeTime
- a piecewise constant function optimized by TreeTime (skyline)

The coalescent contribution to the model is necessary for large trees with many similar sequences.

## Relaxed clock

This feature allows different branches to accumulate mutations at different rates. Each branch has a factor that rescales the mutation rate. These factors are modeled by a probability distribution (either independent on each branch or coupled via an auto-correlated molecular clock) and optimized as part of the iteration loop. For specific probability distributions, the optimal values can be calculated analytically.

## Polytomy resolution

This feature is important for trees with many similar sequences that produce multifurcations. If these multifurcations are resolved into a binary tree randomly (as many tree builders do), the resulting bifurcations will typically be inconsistent with the temporal ordering. TreeTime collapses zero-length branches to create polytomies and then optionally resolves them in a way most consistent with the temporal ordering of the nodes. Two heuristics exist:

- **Greedy**: always merge the pair that increases the likelihood most. This often produces atypical tree shapes.
- **Stochastic**: generates more typical trees while still producing clades compatible with the time ordering.

## Confidence and rate variation

Estimating confidence intervals for internal node dates requires accounting for two sources of noise: the randomness of mutation accumulation, and uncertainty in the rate estimate itself. The forward/backward pass output captures the former (described by the GTR model and branch length distributions). The rate uncertainty is addressed by rerunning the timetree forward/backward pass with rate estimates one standard deviation higher and lower, then combining the central, high, and low estimates.

Uncertainty is output as a table with upper/lower bounds for each node, and can also be specified in the Auspice JSON format.
