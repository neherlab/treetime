# `--gtr-iterations` refits the GTR against frozen marginal messages

`treetime ancestral --model infer --gtr-iterations N` is documented as re-estimating the rate matrix after each reconstruction pass. The code does not re-run reconstruction between fits.

## v1 behavior

- `fn refine_gtr` [packages/treetime/src/gtr/refinement.rs#L17-L88](../../packages/treetime/src/gtr/refinement.rs#L17-L88) with `RefinementRate::Keep` (`GtrRefinement::keep_rate()`) takes the backward and forward messages of the start-model update once, then runs `count_transitions` and `infer_gtr` `N + 1` times (`0..=iterations`) against those same messages, followed by one `marginal_update`
- Each fit changes the GTR through the edge transition matrices, so the iterations are not trivial, but the fixed point depends on the starting messages and is not a stationary point of the tree likelihood. There is no convergence or likelihood check
- The gate is `gtr_iterations > 0 && model == Infer` in [packages/treetime/src/ancestral/plan.rs#L51-L52](../../packages/treetime/src/ancestral/plan.rs#L51-L52), and the refinement runs through `fn fit_marginal()` [packages/treetime/src/ancestral/plan.rs#L292-L309](../../packages/treetime/src/ancestral/plan.rs#L292-L309)

## v0 behavior

- The ancestral CLI fits the GTR once: `infer_gtr` then `_ml_anc` ([packages/legacy/treetime/treetime/treeanc.py#L564-L566](../../packages/legacy/treetime/treetime/treeanc.py#L564-L566)). v0 has no iteration flag
- `infer_gtr_iterative` ([treeanc.py#L1634-L1677](../../packages/legacy/treetime/treetime/treeanc.py#L1634-L1677)) implements a real EM loop with a likelihood stop, but nothing calls it

## Measurement on simulated data

The refined GTR stays biased when the number of sites grows. Simulation: coalescent trees of 40 tips (8 batches of 5 over 2 years, constant `Tc` of 5), a random GTR with `mu = 0.4` per year, exact branch lengths in substitutions per site, so the true rate of the model is 1. Runs of `ancestral::pipeline::run` with `model: GtrModelName::Infer` on the first 1000 and on all 10000 sites:

- Dense, `gtr_iterations: 3`: the rate settles at about 0.3 to 0.4 instead of 1 in the measured replicates, and `pi` moves away from the leaf composition, for example `[0.043, 0.666, 0.170, 0.121]` against the true `[0.058, 0.567, 0.230, 0.146]` at 10000 sites. 20 iterations give almost the same values, so the loop reaches a fixed point that is not the true model
- Dense, `gtr_iterations: 0` (one fit from Fitch counts): the rate is 1 and `pi` is within about 0.02 of the truth
- Paired relation: the log ratio of the error at 10000 and at 1000 sites has a mean of -0.0077 (standard deviation 0.111) for the RMS error of `log w` and -0.0061 (0.158) for the largest error of `pi`, over 200 replicates; one-sided p-values 0.165 and 0.292. 10 times more sites leave both errors almost unchanged
- Sparse, `gtr_iterations: 3` and 20: the rate collapses (about 0.03, then 0.0005) and the exchangeabilities grow without bound; [M-gtr-sparse-transition-counts-diverge-from-dense.md](M-gtr-sparse-transition-counts-diverge-from-dense.md) covers this

Ignored red test: `test_sim_ancestral_gtr_error_falls_with_more_sites` [packages/treetime/src/ancestral/__tests__/test_sim_ancestral_gtr_error_reduction.rs](../../packages/treetime/src/ancestral/__tests__/test_sim_ancestral_gtr_error_reduction.rs) (dense, 200 replicates). The pilots did not show the effect (lower 95% bound of the effect at or below 0 with 200 replicates), so the sizing rule gives no replicate count. After a fix, size the test again from a new pilot and remove the `#[ignore]`.

## Stale documentation

- [kb/decisions/ancestral-iterative-gtr-refinement.md](../decisions/ancestral-iterative-gtr-refinement.md) describes an algorithm (full reconstruction each iteration, with a likelihood stop) that the code does not implement, and no human approval is recorded for it. Its text stays as it is until the question below is decided
- The `--dense` help text ([packages/app-commands/src/commands/ancestral/args.rs#L130-L133](../../packages/app-commands/src/commands/ancestral/args.rs#L130-L133)) and the `--gtr-iterations` help text ([packages/app-commands/src/commands/ancestral/args.rs#L294-L298](../../packages/app-commands/src/commands/ancestral/args.rs#L294-L298)) describe reconstruction passes that do not happen

## Open question

> [!IMPORTANT]
> **Decision required.** The loop is neither v0 behavior nor an EM over the tree likelihood. Options:
>
> - Remove the flag: one GTR fit, as v0 does
> - Replace the loop with a real EM that re-runs the marginal passes each iteration and stops on likelihood change, as the unused v0 `infer_gtr_iterative` does
>
> Evidence: the v1 and v0 behavior above. Mugration uses the same routine with `RefinementRate::Optimize` (v0 order, with a rate step after each fit) and is outside this question.
