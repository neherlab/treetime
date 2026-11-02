# Backward time messages lose accuracy along chains of zero-length branches

Along a chain of internal branches without substitutions, the backward pass of the marginal time inference accumulates a discretization error that grows faster than the chain length. The subtree messages lose mass on the old side, so the node posteriors come out too young and too narrow, and the log marginal likelihood of the time tree is wrong by up to 171 nats on a 200-tip chain. In long chains the posteriors become so narrow that the forward pass stops at the grid point limit ([M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md](M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md)). Tree builders write such chains for sets of identical or nearly identical sequences, and v1 keeps them because it lacks the short-branch collapse of v0 ([M-timetree-pre-loop-short-branch-pruning-missing.md](M-timetree-pre-loop-short-branch-pruning-missing.md)).

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

A sum-product on a fine grid agrees with the closed form to $3 \cdot 10^{-5}$ nats.

## Evidence

With $L = 754$ and $\mu = 0.0014$ ($\lambda = 1.0556$ per year), tips at 2023.0, and the production backward pass:

| Tips | Exact $\log Z$ | v1 $\log Z$ | Difference |
| ---: | ---: | ---: | ---: |
| 20 | -43.364 | -43.248 | 0.116 |
| 50 | -151.129 | -149.236 | 1.893 |
| 200 | -874.000 | -702.484 | 171.516 |

The error of the subtree log mass on the 200-tip chain, against the closed form for the subtree below each node: $2.2 \cdot 10^{-4}$ nats at the lowest node, 0.017 at 10 nodes above it, 1.05 at 40, 35 at 100, and 138 at 180.

On a synthetic 200-tip caterpillar run through the CLI without a coalescent prior, v1 dates `NODE_0000100` in [2022.970, 2023.0], while an exact numeric oracle gives [2022.25, 2022.46], and the root of the backward pass alone is 1.6 years too young. A 400-tip caterpillar stops at the grid point limit with "A grid over [0, 7.20055] with spacing 7.09998e-6 needs 1014165 points"; 50, 100, and 200 tips finish.

Four Nextstrain production builds stop the same way, all with the fixed clock rates of their workflows: rubella E1, oropouche S, and the H5N1 cattle-outbreak builds of PB2 and HA. Their input trees from IQ-TREE hold chains of zero-length internal branches with up to 288 (rubella), 421 (oropouche), 511 (PB2), and 428 (HA) edges; in oropouche, more than 400 of them lie above the earliest tip. Each build finishes at the same fixed rate once the tree is pruned with `treetime prune --collapse-short 0.1/L` before `timetree` (rubella E1 in 14 seconds).

## Mechanism

Near the youngest tip date, the subtree message of a node with $p$ chain nodes below it falls by about $\lambda (p + 1)$ nats per year, so the message gets steeper with every step up the chain. The backward pass keeps the spacing of the branch distribution grid, about 0.023 years in the 200-tip chain, and does not refine it as the message steepens: the convolution of `fn convolve_across_edge()` [packages/treetime/src/timetree/inference/backward_pass.rs#L152](../../packages/treetime/src/timetree/inference/backward_pass.rs#L152) and the products of the child messages discretize a function that changes by several nats within one grid cell. Each step adds an error, and the error per step grows with the slope, so the total grows faster than the chain length.

The forward pass then convolves each too-narrow posterior with a long branch at the posterior's spacing, which requests the grids that exceed the limit (see the narrow-branch issue above).

## Fix direction

The backward messages must resolve their own slope: the spacing of a message must follow the scale over which it changes, as the forward pass refines the spacing of a message whose posterior it under-resolves. The pre-loop collapse of v0 is the reference behavior and removes most such chains from the input trees, but the exact posterior of a chain is well defined, so the backward pass must be correct on any chain that reaches it.

## Reproduction

`packages/treetime/src/timetree/inference/__tests__/test_zero_length_ladder.rs` builds the caterpillar with exact zero-substitution branch likelihoods and runs the production backward and forward passes. Its tests are ignored until the fix:

- `test_zero_length_ladder_log_marginal_likelihood_matches_the_exact_value`, cases `tips_20`, `tips_50`, `tips_200`: compare $\log Z$ with the closed form above within `1e-6` nats
- `test_zero_length_ladder_of_400_tips_dates_every_node_within_the_default_grid_point_limit`: fails at the grid point limit like the production builds; takes about two minutes in the debug build

```bash
./dev/docker/run just test-rs -E 'test(zero_length_ladder)' --run-ignored only
```

## Validation

- The ladder tests pass, and remove their `#[ignore]`
- rubella E1, oropouche S, and the two cattle-outbreak builds finish without the pre-loop collapse at the default grid point limit
- The node posteriors of a 200-tip caterpillar match the exact oracle

## Related

- [M-timetree-time-pass-posteriors-deviate-from-the-exact-posterior.md](M-timetree-time-pass-posteriors-deviate-from-the-exact-posterior.md): the discretization error of the time passes on a three-tip tree, which a chain multiplies
