# Homoplasy Scanner

Identifies recurrent mutations: the same mutation, or mutations at the same site, on several branches of a tree. Recurrence beyond the expectation of independent substitutions points to convergent evolution under selection, recombination, contamination, sequencing artifacts, or sites that evolve fast.

v1: mutation mapping in `fn run_homoplasy()` [`packages/app-commands/src/commands/homoplasy/run.rs`](../../packages/app-commands/src/commands/homoplasy/run.rs), statistics in `fn run()` [`packages/treetime/src/homoplasy/pipeline.rs`](../../packages/treetime/src/homoplasy/pipeline.rs), Poisson comparison in `fn site_histogram()` [`packages/treetime/src/homoplasy/site_hits.rs`](../../packages/treetime/src/homoplasy/site_hits.rs).
v0: `scan_homoplasies()` (`#scan_homoplasies`) in [`packages/legacy/treetime/treetime/wrappers.py#L82-L317`](../../packages/legacy/treetime/treetime/wrappers.py#L82-L317).

Divergences from v0 and their reasons: [kb/decisions/homoplasy-mutation-mapping-and-counting.md](../decisions/homoplasy-mutation-mapping-and-counting.md).

## Algorithm

1. Multiply every input branch length by `--rescale`
2. Reconstruct ancestral sequences with the `ancestral` pipeline (`--method-anc`, marginal by default) and collect the mutations of every branch
3. Sort each mutation with `fn classify_mutation()` [`packages/treetime/src/homoplasy/classify.rs`](../../packages/treetime/src/homoplasy/classify.rs) into substitutions between canonical states, changes involving an ambiguous character, and insertions and deletions. Substitutions are read after unknown-state bridging, so a change through an ancestor with unknown state counts once
4. Group each class by mutation identity and count the branches of each identity; group substitutions also by site
5. Compute the statistics below and rank the identities by multiplicity, then by position

## Statistics of substitutions

Symbols:

- $L$: number of sites, the alignment length plus `--const`
- $M$: number of substitutions, counted once per branch
- $n_k$: number of sites hit by $k$ substitutions, with $n_0 = L - \#\{\text{sites with } k \ge 1\}$
- $K$: largest hit count plus 1
- $t_b$: length of branch $b$

Multiplicity histogram: the number of distinct substitutions $(a, p, d)$ that occur on $m$ branches, for all branches and for terminal branches only.

Tree lengths: the total length is $\sum_b t_b$ over all branches. The terminal length corrects terminal branches for multiple hits with the probability of an odd number of substitutions on a branch,

$$
\sum_{b\ \text{terminal}} e^{-t_b} \sinh t_b = \sum_{b\ \text{terminal}} \tfrac12 \left(1 - e^{-2 t_b}\right)
$$

Poisson comparison: under independent substitutions with a uniform rate, hits per site follow a Poisson distribution with mean $\lambda = M / L$, so

$$
p_k = \frac{\lambda^k e^{-\lambda}}{k!}, \qquad \mathbb{E}[n_k] = L \, p_k
$$

The log-likelihood difference compares the log-likelihood of the observed histogram under this distribution with its expectation under the same distribution:

$$
\Delta = \sum_{k < K} n_k \ln p_k - L \sum_{k < 3K} p_k \ln p_k
$$

A negative $\Delta$ means that substitutions cluster at fewer sites than independent substitutions would produce. With $M = 0$, $n_0 = L$ and $\Delta = 0$.

Per taxon: the substitutions on the terminal branch of each sample at sites hit two or more times anywhere in the tree, so a reversion or another allele at a recurrent site counts.

## Statistics of the other classes

- Ambiguous changes: multiplicity histogram and ranked list by identity, sites with their number of branches, and the number of changes per terminal branch. Many such changes on one terminal branch point to a mixed or low-quality sample
- Insertions and deletions: identity is the kind, the alignment columns, and the inserted or deleted characters. Multiplicity histograms and ranked lists for all and for terminal branches, and the number of terminal indels per sample that occur on two or more branches. No site histogram and no Poisson comparison

## Validation

- Unit tests derive expected counts, Poisson expectations, and $\Delta$ from the formulas above: [`packages/treetime/src/homoplasy/__tests__/test_pipeline.rs`](../../packages/treetime/src/homoplasy/__tests__/test_pipeline.rs)
- Golden master against v0 on `data/zika/86`: [`packages/app-commands/src/commands/homoplasy/__tests__/test_gm_homoplasy.rs`](../../packages/app-commands/src/commands/homoplasy/__tests__/test_gm_homoplasy.rs)
- Dense and sparse reconstruction give equal statistics on `data/zika/86` and `data/flu/h3n2/20`: [`packages/app-commands/src/commands/homoplasy/__tests__/test_dense_sparse.rs`](../../packages/app-commands/src/commands/homoplasy/__tests__/test_dense_sparse.rs). On `data/rsv/a/20` they differ ([kb/issues/M-ancestral-sparse-dense-internal-residues-diverge.md](../issues/M-ancestral-sparse-dense-internal-residues-diverge.md))
