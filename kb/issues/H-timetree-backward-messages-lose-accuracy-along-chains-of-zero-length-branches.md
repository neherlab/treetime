# Backward time messages lose accuracy along chains of zero-length branches

Along a chain of internal branches without substitutions, the backward pass of the marginal time inference accumulates a second-order discretization error that grows faster than the chain length. At the default grid, the log marginal likelihood of the time tree is wrong by 73.9 nats on a 200-tip chain, and the posteriors along such chains are too narrow, which makes the forward pass request fine grids ([M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md](M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md)). Tree builders write such chains for sets of identical or nearly identical sequences. The pre-loop short-branch collapse removes the chains whose sequences agree, as v0 does, but chains remain where the input mode runs without an alignment or the collapse keeps an edge.

## Exact reference

A caterpillar tree of $n$ tips that share one sampling date, with no substitution on any branch, has a closed-form time likelihood. A branch of duration $\tau$ that carries no substitution has the likelihood $e^{-\lambda \tau}$. Order the $n - 1$ internal nodes from the youngest to the root: the interval below each internal node, from the tip date or from the next younger internal node, is crossed by $m$ lineages, with $m$ running from $n$ (below the youngest internal node) down to 2 (below the root). The total branch duration is the sum of the interval lengths weighted by $m$. Without a coalescent prior, the interval lengths are independent, and integrating them out gives

$$\log Z = -\sum_{m=2}^{n} \ln(\lambda m)$$

where:

- $Z$ -- marginal likelihood of the node times, the raw root log mass `log_lh_time` of the time inference
- $\lambda = L \mu$ -- expected substitutions per year on one branch
- $L$ -- sequence length in sites
- $\mu$ -- clock rate in substitutions per site per year
- $m$ -- number of lineages that cross one interval between consecutive internal node times
- $n$ -- number of tips

## Evidence

With $L = 754$ and $\mu = 0.0014$ ($\lambda = 1.0556$ per year), tips at 2023.0, and the production backward pass, the error of $\log Z$ in nats against the number of points of the branch grid:

| Tips | 300 points (default) | 1197 points | 4785 points |
| ---: | ---: | ---: | ---: |
| 20 | 0.116 | 0.0075 | |
| 50 | 1.89 | 0.129 | |
| 200 | 73.9 | 8.13 | 0.546 |

The error falls at second order with the spacing, by 3.75 to 3.96 times per halving, and grows about as $n^3$ at a fixed spacing. The convolution tails below the FFT trust floor are exact up to this discretization error ([kb/decisions/distribution-convolution-tails-by-exponential-tilting.md](../decisions/distribution-convolution-tails-by-exponential-tilting.md)). A 400-tip caterpillar dates every node at the default grid.

Four Nextstrain production builds hold chains of zero-length internal branches in their IQ-TREE input trees, with up to 288 (rubella E1), 421 (oropouche S), 511 (H5N1 cattle-outbreak PB2), and 428 (HA) edges; in oropouche, more than 400 of them lie above the earliest tip. Without the pre-loop collapse, at the fixed clock rates of their workflows, all four stopped at the grid point limit before the convolution tails were computed by exponential tilting; with tilting they finish in 98 (rubella E1), 11 (oropouche S), 58 (HA), and 122 (PB2) seconds in the release build, against 28, 3, 21, and 24 seconds with the collapse. The chains still cost run time, and the margin to the limit is not known.

### Root date against the closed form

The root age of the caterpillar is the second-largest of $n$ independent exponential ages of rate $\lambda$ (the sum of the gaps $m = 2..n$ is that order statistic by the Rényi representation), with the CDF $(1 - u)^{n-1}(1 + (n-1)u)$, $u = e^{-\lambda x}$, the density $\lambda n (n-1) u^2 (1-u)^{n-2}$, the mean $\sum_{m=2}^{n} 1/(\lambda m)$, and the variance $\sum_{m=2}^{n} 1/(\lambda m)^2$. Error of the posterior mean of the root date in years against the closed form, positive when the root is too young; v1 at 300, 1197, and 4785 points of the branch grid (1, 4, and 16 times the default), v0 0.12.1 at its default grid (`treetime --sequence-length 754 --clock-rate 0.0014 --keep-root --keep-polytomies --time-marginal always`, in the Nextstrain image of the pipelandia round):

| Tips | v1 1x | v1 4x | v1 16x | v0 |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 0.0171 | 0.00112 | 0.00007 | 0.046 |
| 50 | 0.106 | 0.0075 | 0.00048 | 0.122 |
| 200 | 0.777 | 0.111 | 0.0078 | 1.287 |

Both engines move the whole posterior toward the tips: the quantiles and the 90% max posterior region move by the same amount as the mean, and the standard deviation is within 0.004 years of the closed form. The v1 error falls at second order once the grid resolves the chain (14 to 16 times per quartering of the spacing). In v0, the rectangle sum of `convolve_fft()` causes the error at 50 tips ([kb/v0-errata/timetree-fft-convolution-rectangle-rule.md](../v0-errata/timetree-fft-convolution-rectangle-rule.md)), but at 200 tips it causes only 0.14 of the 1.29 years.

> [!IMPORTANT]
> **Investigation required.** At 200 tips, 1.15 years of the v0 root date error remain with trapezoid weights in `convolve_fft()`, and they do not fall when $dt$ is divided by 2 or 4 (1.150, 1.189, 1.180 years). No v0 code location for this error is shown, so it is not an erratum. Evidence to collect: the error at 200 tips when each other v0 grid is refined alone, namely the branch grid of `BranchLenInterpolator` (`precision_branch`), the node grids of the backward messages (`Distribution` resampling in `NodeInterpolator.multiply()` and the `effective_support` that crops the operands of `convolve_fft()`), and the cut at `1e-13` of the peak with its straight-line tails.

## Mechanism

Near the youngest tip date, the subtree message of a node with $p$ chain nodes below it falls by about $\lambda (p + 1)$ nats per year, so the message gets steeper with every step up the chain. The backward pass keeps the spacing of the branch distribution grid, about 0.023 years in the 200-tip chain, and does not refine it as the message steepens: the convolution of `fn convolve_across_edge()` [packages/treetime/src/timetree/inference/backward_pass.rs#L152](../../packages/treetime/src/timetree/inference/backward_pass.rs#L152) and the products of the child messages discretize a function that changes by several nats within one grid cell. Each step adds a second-order error, and the error per step grows with the slope, so the total grows faster than the chain length.

The forward pass then convolves each too-narrow posterior with a long branch at the posterior's spacing, which requests the grids that exceed the limit (see the narrow-branch issue above).

## Fix direction

The backward messages must resolve their own slope: the spacing of a message must follow the scale over which it changes, as the forward pass refines the spacing of a message whose posterior it under-resolves. A prototype that refined the spacing of each backward message to its mass width reduced the ladder errors by 6 and 4 times, but requested 1.4 million points on oropouche S, so the spacing rule must also bound the grid. The pre-loop collapse of v0 is the reference behavior and removes most such chains from the input trees, but the exact posterior of a chain is well defined, so the backward pass must be correct on any chain that reaches it.

## Reproduction

`packages/treetime/src/timetree/inference/__tests__/test_zero_length_ladder.rs` builds the caterpillar with exact zero-substitution branch likelihoods and runs the production backward and forward passes:

- `test_zero_length_ladder_log_marginal_likelihood_converges_at_second_order_to_the_exact_value`: the error of $\log Z$ against the closed form falls by at least 3 times per halving of the branch grid spacing
- `test_zero_length_ladder_root_date_posterior_converges_to_the_hypoexponential_root_age`: the largest error of the mean, standard deviation, 5% and 95% quantiles, and 90% max posterior region of the root date against the closed form falls with each quartering of the spacing, by at least 9 times on the last
- `test_zero_length_ladder_of_400_tips_dates_every_node_within_the_default_grid_point_limit`

## Validation

- The error of $\log Z$ on the 200-tip caterpillar at the default grid falls well below its value of 73.9 nats
- rubella E1, oropouche S, and the two cattle-outbreak builds still finish without the pre-loop collapse at the default grid point limit

## Related

- [M-timetree-time-pass-posteriors-deviate-from-the-exact-posterior.md](M-timetree-time-pass-posteriors-deviate-from-the-exact-posterior.md): the discretization error of the time passes on a three-tip tree, which a chain multiplies
- [N-timetree-resampled-time-distributions-with-two-hard-ends-move-one-end-by-up-to-half-a-cell.md](N-timetree-resampled-time-distributions-with-two-hard-ends-move-one-end-by-up-to-half-a-cell.md): a resampled operand with two hard ends, such as the date range of a tip dated to a year, moves one hard end by up to half a cell
