# Iterative outer GTR refinement for ancestral reconstruction

This document proposes an optional extension to the ancestral command: when running marginal ancestral reconstruction with `--model infer`, alternate between ancestral-state reconstruction and GTR re-estimation for several outer iterations instead of performing only one conditional model fit.

The proposal is **not accepted** and **not implemented**. It records the scientific case, the relevant current code, and the main tradeoffs.

It also records the current accepted behavior that motivated the proposal: both v0 and v1 ancestral command paths use one-shot conditional GTR inference, while v0 mugration uses a true outer refinement loop.

## Current behavior

### v0 ancestral CLI

The v0 ancestral wrapper at [`packages/legacy/treetime/treetime/wrappers.py#L588-L623`](../../packages/legacy/treetime/treetime/wrappers.py#L588-L623) calls:

```python
treeanc.infer_ancestral_sequences(
    'ml',
    infer_gtr=params.gtr == 'infer',
    marginal=params.marginal,
    fixed_pi=fixed_pi,
    reconstruct_tip_states=params.reconstruct_tip_states,
)
```

Inside [`TreeAnc.infer_ancestral_sequences()`](../../packages/legacy/treetime/treetime/treeanc.py#L516-L565), `infer_gtr=True` triggers one call to [`TreeAnc.infer_gtr()`](../../packages/legacy/treetime/treetime/treeanc.py#L1500-L1589) before the actual ancestral reconstruction pass.

The ancestral CLI does not call [`TreeAnc.infer_gtr_iterative()`](../../packages/legacy/treetime/treetime/treeanc.py#L1634-L1677) and does not call [`TreeAnc.optimize_gtr_rate()`](../../packages/legacy/treetime/treetime/treeanc.py#L1679-L1704).

### v1 ancestral CLI

The v1 ancestral command in [`packages/treetime/src/commands/ancestral/run.rs`](../../packages/treetime/src/commands/ancestral/run.rs) also performs one conditional model fit.

- Sparse mode: compress sequences, infer one model with [`get_gtr_sparse()`](../../packages/treetime/src/gtr/get_gtr.rs#L115-L132), then run one marginal reconstruction under that model.
- Dense mode: run one preliminary marginal pass under temporary JC69 to populate profiles, infer one model with [`get_gtr_dense()`](../../packages/treetime/src/gtr/get_gtr.rs#L135-L152), then rerun marginal under the inferred model.

Relevant orchestration is at [`packages/treetime/src/commands/ancestral/run.rs#L109-L188`](../../packages/treetime/src/commands/ancestral/run.rs#L109-L188).

This one-shot behavior is the current accepted ancestral workflow. The proposal does not claim that current v1 is missing an already-accepted algorithm. It proposes a possible future extension beyond the currently accepted workflow.

### Mugration contrast

The mugration wrapper in v0 does implement an outer refinement loop. [`reconstruct_discrete_traits()`](../../packages/legacy/treetime/treetime/wrappers.py#L653-L811) performs:

1. initial `infer_ancestral_sequences(... infer_gtr=True, marginal=True)`
2. `optimize_gtr_rate()`
3. repeated `infer_gtr(...)` plus `optimize_gtr_rate()`
4. final `infer_ancestral_sequences(infer_gtr=False, marginal=True)`

This difference is the direct motivation for the proposal.

## Scientific background

### Substitution models on trees

Ancestral sequence reconstruction under GTR treats sequence evolution along each branch as a continuous-time Markov chain with transition matrix

$$
P(t) = e^{Qt},
$$

where `Q` is the instantaneous rate matrix. In TreeTime's GTR parameterization, `Q` is determined by:

- equilibrium frequencies `pi`
- symmetric exchangeability matrix `W`
- overall rate scalar `mu`

This is the standard finite-state phylogenetic model family described by Tavaré (1986) and used by Felsenstein-style likelihood calculations.

### Why outer refinement is scientifically plausible

Internal-node sequences are latent variables. When `--model infer` is active, the ancestral command is trying to optimize a likelihood over both:

- hidden ancestral states or their posterior distributions
- model parameters `W`, `pi`, and `mu`

That is a latent-variable maximum-likelihood problem. Alternating between:

1. reconstructing or marginalizing ancestral states under the current model, and
2. re-estimating model parameters from the resulting expected transition counts,

is an EM-like or coordinate-ascent strategy. It is scientifically standard to expect such alternation to improve or stabilize the fit when the current model is far from the local optimum.

The current v1 GTR core already performs **inner** iteration on `W`, `pi`, and `mu` for fixed sufficient statistics in [`infer_gtr_impl()`](../../packages/treetime/src/gtr/infer_gtr/common.rs#L97-L157). The proposal adds **outer** iteration over the latent-state statistics themselves.

### Inner iteration vs outer refinement

The distinction between the accepted current workflow and the proposed extension is central.

- Inner GTR iteration solves for `W`, `pi`, and `mu` while holding expected mutation counts and time-in-state totals fixed.
- Outer refinement recomputes those sufficient statistics after ancestral posteriors change, then re-estimates the model, then reconstructs again.

In other words:

1. current ancestral workflow = one outer conditional update plus an inner parameter solver
2. proposed workflow = repeated outer conditional updates, each of which contains the same inner parameter solver

That difference is why the material belongs in a proposal rather than in the accepted algorithm inventory.

### Why the expected benefit is smaller than for mugration

Mugration effectively reconstructs one discrete character over the tree. In that setting, small shifts in `pi`, `W`, or `mu` can change the posterior at ambiguous internal nodes substantially.

Sequence ancestral reconstruction aggregates information over many sites. The expected substitution counts and time-in-state totals are therefore usually better constrained. A single conditional GTR fit often captures most of the available signal.

The proposal is still scientifically plausible, but the expected payoff is concentrated in cases such as:

- short alignments
- strong compositional bias
- large alphabets or amino-acid models
- weakly resolved internal branches
- datasets where the initialization model is far from the fitted optimum

## Related current code

### v1 model inference core

- [`packages/treetime/src/gtr/infer_gtr/common.rs#L97-L157`](../../packages/treetime/src/gtr/infer_gtr/common.rs#L97-L157): iterative update of `W`, `pi`, and `mu` for fixed counts
- [`packages/treetime/src/gtr/infer_gtr/dense.rs#L15-L23`](../../packages/treetime/src/gtr/infer_gtr/dense.rs#L15-L23): dense inference entry point
- [`packages/treetime/src/gtr/infer_gtr/sparse.rs#L10-L18`](../../packages/treetime/src/gtr/infer_gtr/sparse.rs#L10-L18): sparse inference entry point

### v1 marginal passes

- [`packages/treetime/src/commands/ancestral/marginal.rs#L17-L49`](../../packages/treetime/src/commands/ancestral/marginal.rs#L17-L49): `initialize_marginal()` and `update_marginal()`
- [`packages/treetime/src/representation/partition/marginal_dense.rs`](../../packages/treetime/src/representation/partition/marginal_dense.rs)
- [`packages/treetime/src/representation/partition/marginal_passes.rs`](../../packages/treetime/src/representation/partition/marginal_passes.rs)

### v0 precedent

- [`packages/legacy/treetime/treetime/treeanc.py#L1634-L1677`](../../packages/legacy/treetime/treetime/treeanc.py#L1634-L1677): `infer_gtr_iterative()` exists as a reusable method
- [`packages/legacy/treetime/treetime/wrappers.py#L653-L811`](../../packages/legacy/treetime/treetime/wrappers.py#L653-L811): mugration uses outer refinement in production CLI code

## Proposal

### Scope

Add optional outer refinement for:

- `ancestral`
- `--method-anc marginal`
- `--model infer`

Do not extend parsimony or joint reconstruction in this proposal.

### Proposed algorithm

For either dense or sparse marginal ancestral reconstruction:

1. Initialize with the current startup model, matching existing behavior.
2. Run marginal reconstruction to obtain current posteriors and log likelihood.
3. Infer a new GTR from the current sufficient statistics.
4. If the new model changes materially or the log likelihood still improves, rerun marginal reconstruction under the new model.
5. Stop when improvement falls below a tolerance or after a small fixed number of iterations.

This would convert the current one-shot conditional update into a bounded outer alternation.

### Recommended first version

The first version should **not** mirror mugration mechanically.

Recommended differences from mugration:

- reuse v1's existing `infer_gtr_impl()` for parameter updates
- add an outer loop around `update_marginal()` plus `get_gtr_dense()` or `get_gtr_sparse()`
- use likelihood-based stopping criteria
- avoid adding a separate Brent-style `mu` optimization step unless validation demonstrates a measurable benefit

This recommendation follows from the scientific role of the data:

- mugration has one effective site, so the overall switching rate is weakly identified and extra `mu` optimization is more likely to matter
- sequence ancestral reconstruction already estimates `mu` from large multi-site counts inside `infer_gtr_impl()`

## Expected benefits

### Potential gains

- improved fit when the initial model is far from the local optimum
- better `pi` estimates on compositionally biased datasets
- fewer model-dependent posterior artifacts at ambiguous internal nodes
- closer alignment between latent-state inference and model estimation

### Why this is scientifically defensible

The proposal aligns the optimization target more closely with the actual latent-variable likelihood. It reduces the mismatch between:

- counts computed under one model, and
- final reconstructed states produced under another model.

One-shot inference leaves that mismatch in place after the final rerun. Outer refinement reduces it by recomputing counts after the model changes.

## Risks and caveats

### Overfitting and numerical churn

Extra outer iterations can chase small likelihood changes that do not matter biologically. On typical multi-site nucleotide datasets, the improvement might be negligible while increasing runtime and numerical sensitivity.

### Posterior-count feedback

If posterior profiles are overconfident, repeated re-estimation can reinforce early biases. This is a standard latent-variable optimization risk. It argues for:

- marginal, not joint, updates
- bounded iteration count
- explicit stopping rules
- parity and stability validation

### v0 parity

This proposal would not necessarily make v1 more like the v0 ancestral CLI, because the v0 ancestral CLI also stops after one conditional fit. It is therefore best framed as a **candidate intentional extension**, not as a parity task.

## Current conclusion motivating the proposal

The codebase currently supports the following scientifically meaningful distinction:

- v0 ancestral: one-shot conditional GTR inference
- v1 ancestral: one-shot conditional GTR inference
- v0 mugration: iterative outer GTR refinement

This proposal exists because mugration demonstrates that outer refinement is scientifically plausible in the project domain, while ancestral currently uses the simpler one-shot workflow.

## Validation plan if accepted

1. Add example-based tests showing monotonic or non-decreasing sequence log likelihood across outer iterations, within numerical tolerance.
2. Compare one-shot and iterative outputs on datasets with strong compositional skew.
3. Check whether inferred `pi`, `W`, and `mu` stabilize within a small number of iterations.
4. Confirm that runtime remains acceptable on representative large datasets.
5. Decide whether the feature is always-on or exposed behind a dedicated CLI flag.

## Recommendation

This proposal is scientifically plausible and internally consistent with TreeTime's existing GTR machinery. It is worth considering as an intentional extension for `ancestral --model infer`, but it should be treated as an optimization-quality improvement rather than a parity requirement.

The main reason to accept it is improved joint consistency between latent ancestral profiles and inferred substitution parameters. The main reason to reject it is that the expected practical benefit for long, informative alignments may be too small to justify the extra complexity.

## References

- Felsenstein, J. (1981). "Evolutionary trees from DNA sequences: A maximum likelihood approach." _Journal of Molecular Evolution_, 17(6), 368-376.
- Tavaré, S. (1986). "Some probabilistic and statistical problems in the analysis of DNA sequences." _Lectures on Mathematics in the Life Sciences_, 17, 57-86.
- Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). "Maximum likelihood from incomplete data via the EM algorithm." _Journal of the Royal Statistical Society: Series B_, 39(1), 1-38.
