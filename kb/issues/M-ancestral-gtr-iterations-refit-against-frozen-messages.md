# `--gtr-iterations` refits the GTR against frozen marginal messages

`treetime ancestral --model infer --gtr-iterations N` is documented as re-estimating the rate matrix after each reconstruction pass. The code does not re-run reconstruction between fits.

## v1 behavior

- `fn refine_gtr_model` [packages/treetime/src/gtr/refinement.rs#L18-L48](../../packages/treetime/src/gtr/refinement.rs#L18-L48) destructures the backward and forward messages once, then runs `count_transitions` and `infer_gtr` `N + 1` times (`0..=iterations`) against those same messages, followed by one `marginal_update`
- Each fit changes the GTR through the edge transition matrices, so the iterations are not trivial, but the fixed point depends on the starting messages and is not a stationary point of the tree likelihood. There is no convergence or likelihood check
- The gate is `gtr_iterations > 0 && model == Infer` in [packages/treetime/src/ancestral/plan.rs#L69-L70](../../packages/treetime/src/ancestral/plan.rs#L69-L70)

## v0 behavior

- The ancestral CLI fits the GTR once: `infer_gtr` then `_ml_anc` ([packages/legacy/treetime/treetime/treeanc.py#L564-L566](../../packages/legacy/treetime/treetime/treeanc.py#L564-L566)). v0 has no iteration flag
- `infer_gtr_iterative` ([treeanc.py#L1634-L1677](../../packages/legacy/treetime/treetime/treeanc.py#L1634-L1677)) implements a real EM loop with a likelihood stop, but nothing calls it

## Stale documentation

- [kb/decisions/ancestral-iterative-gtr-refinement.md](../decisions/ancestral-iterative-gtr-refinement.md) describes functions that no longer exist and an algorithm (full reconstruction each iteration) that the code does not implement; no human approval is recorded for it
- The `--gtr-iterations` and `--dense` help texts in [packages/app-commands/src/commands/ancestral/args.rs](../../packages/app-commands/src/commands/ancestral/args.rs) describe reconstruction passes that do not happen

## Open question

Remove the flag (v0 parity, one fit), or replace the loop with a real EM that re-runs the marginal passes each iteration and stops on likelihood change.
