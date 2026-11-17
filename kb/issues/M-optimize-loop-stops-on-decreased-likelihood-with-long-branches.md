# The optimize loop often stops on a decreased likelihood with long branches

The optimize loop updates the branch lengths of all edges from the messages of the previous iteration, and stops when the log-likelihood decreases, returning the iteration with the highest likelihood ([kb/decisions/optimize-damping-disabled-by-default.md](../decisions/optimize-damping-disabled-by-default.md)). On the example datasets every run stopped on convergence or after a single decline.

On simulated JC69 data with continuous gamma site rates (α = 0.5, 10000 sites) on the tree of `data/zika/86` with branch lengths times 10 (0.66 substitutions per site), optimize with `--model jc69` and without topology steps stopped on a decreased log-likelihood in:

- 4 of 20 replicates without rate variation, and 10 more stopped as oscillating
- 14 of 20 replicates with `--partition-by rate-categories`
- 14 of 20 replicates with `--gamma-categories 4`

On the unscaled tree no replicate stopped this way. The returned branch lengths belong to the best iteration, so the result is valid, but the loop may stop before it converges, and the stop rate grows with rate variation.

> [!IMPORTANT]
> **Investigation required.** Record the log-likelihood of each iteration of these runs and the size of the decrease. Check whether `--damping` removes the stops, and whether the decrease comes from the simultaneous update of all edges or from the sparse approximation (repeat in dense mode). The bias measurement of the rate buckets in [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md#partition-schemes) uses these runs.

## Validation

- On the scaled tree, optimize stops on convergence in at least 19 of 20 replicates in each of the three settings, or the decision on damping records why a stop on a decreased likelihood is the intended result
