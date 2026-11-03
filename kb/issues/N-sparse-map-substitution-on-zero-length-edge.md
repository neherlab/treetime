# Sparse marginal reconstruction places substitutions on edges of length 0

At branch length 0 the transition matrix is the identity, so the exact marginal posteriors of a parent and its child are equal, and the most likely states of the two nodes agree at every site. The sparse backend still reports substitutions on some edges of length 0. `fn compute_ml_subs_for_nodes()` [packages/treetime/src/partition/marginal/sparse/forward.rs#L334-L364](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L334-L364) reads the most likely state of each node separately (`fn map_state()`) and reports every site where the two differ.

Observed with sparse `optimize` on `rust` at `a8f16336c`, before the sparse missing-data fix:

| Input | Edges of length 0 with a marginal substitution |
|---|---|
| ncov build input of the Nextstrain workflow (`--model infer --no-indels`) | four edges, with `T28093C`, `T2470C`, `T26936C`, `T23075C` |
| `flu/h3n2/500`, 7 siblings that share `T906C` merged by hand under a new node | the new edge, which optimizes to length 0, carries `T906C`; dense mode places no substitution there |

The topology cleanup of `optimize` reads these substitutions: the short-edge collapse after the loop keeps an internal edge that changes a state (`fn collapse_short_internal_edges()` [packages/treetime/src/optimize/topology/collapse.rs#L84](../../packages/treetime/src/optimize/topology/collapse.rs#L84)), and merges and hoists count a Fitch substitution only where the marginal reconstruction also places it (`fn SparseReconstruction::supported_fitch_subs()` [packages/treetime/src/partition/marginal/reconstruction.rs#L352](../../packages/treetime/src/partition/marginal/reconstruction.rs#L352)). Similar edges were seen on the `dengue`, `flu` and `rsv` datasets.

> [!IMPORTANT]
> **Investigation required.** The observations come from code before the sparse missing-data fix and are not checked again. Rerun sparse and dense `optimize` on `flu/h3n2/500` and on an input with long `N` runs, list the internal edges of length 0 with a marginal substitution, and compare the posteriors of both end nodes at each such site. Candidate causes: the folding of resolved sites into the per-state profiles ([kb/decisions/ancestral-marginal-sparse-parsimony-chain.md](../decisions/ancestral-marginal-sparse-parsimony-chain.md), posterior above 1 - 1e-4), which can fold the child's site while the parent keeps a per-site profile, and the sparse and dense difference of [M-sparse-dense-likelihood-residual-without-missing-data.md](M-sparse-dense-likelihood-residual-without-missing-data.md).
