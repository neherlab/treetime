# Full forward-backward reconstruction per mugration GTR iteration

## Background

GTR parameter estimation from phylogenetic data follows the Expectation-Maximization (EM) framework (Dempster, Laird & Rubin, 1977). EM alternates between computing posterior state distributions given current parameters (E-step) and re-estimating parameters from those posteriors (M-step). For GTR on a tree, the E-step is Felsenstein's pruning algorithm (forward-backward message passing), and the M-step estimates substitution rates and equilibrium frequencies from expected transition counts (Holmes & Rubin, 2002).

Convergence guarantees require that the E-step uses messages computed under a single consistent model. Mixing messages from different models (stale forward + fresh backward) breaks the EM monotonicity property: log-likelihood is no longer guaranteed to increase each iteration.

## Current behavior

The iterative GTR refinement in `refine_gtr_iterative()` ([packages/treetime/src/commands/mugration/gtr_refinement.rs](../../packages/treetime/src/commands/mugration/gtr_refinement.rs)) uses stale forward-pass messages (`msg_to_child`) from the initial reconstruction. Each iteration updates backward messages via rate optimization but does not refresh forward messages. This matches v0's behavior: `infer_gtr(marginal=True)` ([packages/legacy/treetime/treetime/treeanc.py#L1543-L1545](../../packages/legacy/treetime/treetime/treeanc.py#L1543-L1545)) skips `_ml_anc_marginal()` when `sequence_reconstruction == 'marginal'`.

The mutation counting formula in `get_branch_mutation_matrix()` ([packages/treetime/src/gtr/infer_gtr/dense.rs#L63](../../packages/treetime/src/gtr/infer_gtr/dense.rs#L63)) computes the joint `P(child=i, parent=j | data)` from three components that may be under different GTR models:

- `msg_to_parent[i]`: subtree likelihood from backward pass with current GTR (fresh)
- `exp_qt[i,j]`: transition matrix from current GTR (fresh)
- `msg_to_child[j]`: outgroup likelihood from forward pass with initial GTR (stale)

## Proposed change

Insert `run_discrete_marginal(graph, partition)?` before each `count_transitions_discrete()` in the iterative loop:

```rust
for i in 0..iterations {
    run_discrete_marginal(graph, partition)?;
    let counts = count_transitions_discrete(graph, partition)?;
    // ...
}
```

## Expected impact

- Correct EM iterations: both forward and backward messages computed under the current model
- More iterations needed for convergence (each iteration starts fresh)
- Different fixed point from both current v1 and v0 (neither runs correct EM)
- All golden master tests would need recapture (v0 oracles are no longer the target)
- Computational cost: ~2x per iteration (full forward-backward vs backward-only)

## Validation plan

1. Compare v1 (current, stale forward) vs v1 (proposed, fresh forward) on all 7 datasets
2. Check convergence: log-likelihood should monotonically increase across iterations
3. Check stability: increasing iterations beyond 5 should not change results
4. Independent validation against BEAST or PhyML mugration output on a test dataset

## References

- Dempster, A.P., Laird, N.M. & Rubin, D.B. (1977). Maximum likelihood from incomplete data via the EM algorithm. _J R Stat Soc B_, 39(1):1-38. [doi:10.1111/j.2517-6161.1977.tb01600.x](https://doi.org/10.1111/j.2517-6161.1977.tb01600.x)
- Holmes, I. & Rubin, G.M. (2002). An expectation maximization algorithm for training hidden substitution models. _J Mol Biol_, 317(5):753-764. [PMID:11955022](https://pubmed.ncbi.nlm.nih.gov/11955022/)
- Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. _J Mol Evol_, 17(6):368-376. [PMID:7288891](https://pubmed.ncbi.nlm.nih.gov/7288891/)
