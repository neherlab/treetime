# Mugration runtime grows with the cube of the state count

`mugration` with a trait of a few hundred states takes many minutes, because every marginal pass forms the full transition matrix of every edge.

## Measurements

Release build (`just build-release treetime`), `sc2/4500` (4500 tips), default settings, on a 20-thread i9-10900K that other jobs shared (load average about 20 at the start of each run):

| Attribute        | States | Wall time | CPU time | Peak memory |
| ---------------- | ------ | --------- | -------- | ----------- |
| `country`        | 95     | 25 s      | 227 s    | 206 MB      |
| `division`       | 261    | 501 s     | 6769 s   | 535 MB      |
| `submitting_lab` | 321    | 1097 s    | 15008 s  | 551 MB      |

From 95 to 261 states (2.75 times as many) the wall time grows 20 times, which is the cube of 2.75. A trait with about 700 states, such as `pango_lineage` of the same dataset, would take about 3 hours at this rate.

## Cause

- The backward and forward passes call `GTR::propagate_profile()` and `GTR::evolve()` for each edge ([packages/treetime/src/partition/marginal/shared/pass.rs#L196](../../packages/treetime/src/partition/marginal/shared/pass.rs#L196), [#L327](../../packages/treetime/src/partition/marginal/shared/pass.rs#L327)). Without site rates, both form the full matrix `exp(Qt)` through `fn expQt_with_rate()` ([packages/treetime/src/gtr/gtr.rs#L147-L152](../../packages/treetime/src/gtr/gtr.rs#L147-L152)), a product of two n by n matrices, so each edge costs O(n³) for n states
- A discrete trait has one site, so the message of an edge is one row of n values. The product of that row with the matrix costs only O(n²)
- The rate step of the GTR refinement evaluates its cost function about 20 to 30 times per iteration, each time with a full backward pass, so the backward passes dominate the run
- The transition counts of each refinement iteration need every entry of `exp(Qt)` ([packages/treetime/src/partition/marginal/shared/data.rs#L36](../../packages/treetime/src/partition/marginal/shared/data.rs#L36)), so that step stays O(n³) per edge, but it runs once per iteration only

v0 has the same order of cost: its `propagate_profile()` and `evolve()` also form `expQt(t)`, with two n by n products (`np.diag`, then `v.dot(eLambdaT.dot(v_inv))`, [packages/legacy/treetime/treetime/gtr.py#L965-L1066](../../packages/legacy/treetime/treetime/gtr.py#L965-L1066)).

## Fix direction

Propagate messages with fewer rows than states in the eigenbasis: multiply the profile with the eigenvectors, scale each column by `exp(λt)`, and multiply with the inverse eigenvectors (`profile.dot(&v)`, the scaling, then `.dot(&v_inv)` for `propagate_profile()`), which costs O(n²) per row. The site-rate branch of the same functions already computes this way ([packages/treetime/src/gtr/gtr.rs#L160-L164](../../packages/treetime/src/gtr/gtr.rs#L160-L164)). Dense nucleotide alignments (many rows, 4 states) keep the matrix product.

> [!IMPORTANT]
> **Investigation required.** The eigenbasis product clamps negative values of the result, while `fn expQt_with_rate()` clamps negative entries of the matrix, and the two paths round differently. Measure how far the discrete outputs move, and check the mugration golden masters (tolerance 1e-7 against v0), before the change is adopted.

## Validation

- The three runs above, timed again on the same host
- Mugration golden masters within 1e-7 of v0
- The posteriors and log likelihood of the discrete partition agree with the current matrix path within a few ulps on small trees with 2 to 300 states
