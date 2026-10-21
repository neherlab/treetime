# `--gtr-iterations` refits the GTR against frozen marginal messages

`treetime ancestral --model infer --gtr-iterations N` is documented as re-estimating the rate matrix after each reconstruction pass. The code does not re-run reconstruction between fits.

## v1 behavior

- `fn refine_gtr` [packages/treetime/src/gtr/refinement.rs#L17-L88](../../packages/treetime/src/gtr/refinement.rs#L17-L88) with `RefinementRate::Keep` (`GtrRefinement::keep_rate()`) takes the backward and forward messages of the start-model update once, then runs `count_transitions` and `infer_gtr` `N + 1` times (`0..=iterations`) against those same messages, followed by one `marginal_update`
- Each fit changes the GTR through the edge transition matrices, so the iterations are not trivial, but the fixed point depends on the starting messages and is not a stationary point of the tree likelihood. There is no convergence or likelihood check
- The gate is `gtr_iterations > 0 && model == Infer` in [packages/treetime/src/ancestral/plan.rs#L51-L52](../../packages/treetime/src/ancestral/plan.rs#L51-L52), and the refinement runs through `fn fit_marginal()` [packages/treetime/src/ancestral/plan.rs#L292-L309](../../packages/treetime/src/ancestral/plan.rs#L292-L309)

## v0 behavior

- The ancestral CLI fits the GTR once: `infer_gtr` then `_ml_anc` ([packages/legacy/treetime/treetime/treeanc.py#L564-L566](../../packages/legacy/treetime/treetime/treeanc.py#L564-L566)). v0 has no iteration flag
- `infer_gtr_iterative` ([treeanc.py#L1634-L1677](../../packages/legacy/treetime/treetime/treeanc.py#L1634-L1677)) implements a real EM loop with a likelihood stop, but nothing calls it

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
