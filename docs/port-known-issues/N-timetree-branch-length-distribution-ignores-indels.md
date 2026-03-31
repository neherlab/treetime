# Timetree branch length distribution ignores indels

The timetree command's branch length distribution calculation at [packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L50](../../packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L50) calls `evaluate_mixed_log_lh_only()` directly, which evaluates substitution likelihood only. The Poisson indel contribution added to `run_optimize_mixed()` does not reach this code path.

## Current behavior

`compute_branch_length_distribution()` builds a grid of branch lengths and evaluates the substitution log-likelihood at each point via `evaluate_mixed_log_lh_only()`. This produces a branch length distribution used for time inference. Indel events on the edge do not influence this distribution.

## Expected behavior

The branch length distribution should include the indel contribution, consistent with how `run_optimize_mixed()` evaluates edges. An edge with indels but no substitutions should have a distribution peaked away from zero.

## Impact

Low for typical viral datasets (indels are rare). Visible on branches where the only evolutionary signal is an indel event: the timetree distribution will peak at zero while the optimize command correctly assigns positive length.

## Fix

Pass the indel count and rate to `compute_branch_length_distribution()` and use `evaluate_with_indels_log_lh_only()` (or equivalent) instead of `evaluate_mixed_log_lh_only()`.

## Broader design question

The Poisson indel model (documented in [docs/algorithms/optimize.md#L31-L46](../../docs/algorithms/optimize.md#L31-L46)) currently contributes to:

1. Newton-step branch length optimization in `run_optimize_mixed()` -- implemented
2. Branch length distribution grid in timetree -- **not implemented** (this issue)
3. Full branch length likelihood evaluation -- **not implemented**

The `2026-03-24` design note raises extending indel contributions to the full likelihood: "An outstanding question is whether we can extend this to indels. [...] A branch with no subs but an indel should have a finite length. But for this to happen we need to add it somehow to the LH."

The Poisson model ($\ell_{\text{indel}}(t) = k \ln(\mu t) - \mu t - \ln(k!)$) treats each indel event with equal weight regardless of indel length and does not distinguish insertions from deletions. More sophisticated probabilistic indel models exist <a id="cite-1"></a>[Thorne, Kishino, and Felsenstein 1991](https://doi.org/10.1007/BF02193625) [[1](#ref-1)]; <a id="cite-2"></a>[Thorne, Kishino, and Felsenstein 1992](https://doi.org/10.1007/BF00163848) [[2](#ref-2)]; <a id="cite-3"></a>[Bouchard-Cote and Jordan 2013](https://doi.org/10.1073/pnas.1220450110) [[3](#ref-3)]; <a id="cite-4"></a>[Miklos, Lunter, and Holmes 2004](https://doi.org/10.1093/molbev/msh043) [[4](#ref-4)] but add substantial complexity. The current Poisson model is appropriate for the branch-length-prevents-zero use case. Whether it belongs in the full likelihood depends on whether indel events should influence the overall tree likelihood or only edge-level optimization.

## References

1. <a id="ref-1"></a> Thorne, Jeffrey L., Hirohisa Kishino, and Joseph Felsenstein. 1991. "An Evolutionary Model for Maximum Likelihood Alignment of DNA Sequences." _Journal of Molecular Evolution_ 33(2):114-124. https://doi.org/10.1007/BF02193625 [↩](#cite-1)
2. <a id="ref-2"></a> Thorne, Jeffrey L., Hirohisa Kishino, and Joseph Felsenstein. 1992. "Inching toward Reality: An Improved Likelihood Model of Sequence Evolution." _Journal of Molecular Evolution_ 34(1):3-16. https://doi.org/10.1007/BF00163848 [↩](#cite-2)
3. <a id="ref-3"></a> Bouchard-Cote, Alexandre, and Michael I. Jordan. 2013. "Evolutionary Inference via the Poisson Indel Process." _Proceedings of the National Academy of Sciences_ 110(4):1160-1166. https://doi.org/10.1073/pnas.1220450110 [↩](#cite-3)
4. <a id="ref-4"></a> Miklos, Istvan, Gerton A. Lunter, and Ian Holmes. 2004. "A 'Long Indel' Model for Evolutionary Sequence Alignment." _Molecular Biology and Evolution_ 21(3):529-540. https://doi.org/10.1093/molbev/msh043 [↩](#cite-4)
