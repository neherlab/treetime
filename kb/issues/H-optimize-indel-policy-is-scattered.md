# Indel inclusion policy is scattered through optimization and timetree layers

`no_indels: bool` is carried from command arguments through optimize and timetree parameters, initial guesses, likelihood composition, edge optimization, confidence estimation, and branch-distribution inference [`packages/treetime/src/optimize/pipeline.rs#L47`](../../packages/treetime/src/optimize/pipeline.rs#L47) [`packages/treetime/src/optimize/run_loop.rs#L70`](../../packages/treetime/src/optimize/run_loop.rs#L70) [`packages/treetime/src/timetree/inference/runner.rs#L74`](../../packages/treetime/src/timetree/inference/runner.rs#L74) [`packages/treetime/src/timetree/confidence.rs#L59`](../../packages/treetime/src/timetree/confidence.rs#L59).

Receivers independently interpret the flag to control rate estimation, likelihood contribution, zero-branch validity, initial guesses, and distributions. A change to indel treatment therefore requires coordinated semantic edits across unrelated layers.

## Control path

The value flows from command arguments into `OptimizeParams` or `TimetreeParams`, then through optimization loops, dispatch, timetree refinement, confidence estimation, and inference. Its meaning is reconstructed at each layer rather than carried by one policy object.

The affected behaviors are coupled scientifically:

- whether indel rates are estimated;
- whether Poisson indel likelihood contributes to the objective;
- whether zero branch length is admissible;
- which initial guess and branch distribution are constructed;
- which confidence calculation is reported.

Parse command input into one typed indel-likelihood policy that owns these behaviors. The approved Poisson model and the existing include/exclude semantics remain unchanged.

## Required boundary

The policy should expose domain operations or construct the evaluator needed by downstream algorithms. Receivers must not branch independently on a copied boolean. If optimize and timetree require different projections of the policy, derive both from the same parsed state.

## Validation

- Golden-master comparisons cover include/exclude modes in optimize and timetree.
- Unit tests verify rate estimation, objective composition, zero-boundary behavior, initial guesses, and distributions for each policy variant.
- Semantic search finds no forwarded `no_indels: bool` below the parsing boundary.
- Error messages and CLI behavior remain unchanged unless separately approved.

## Related decisions

- [kb/decisions/optimize-indel-contribution-to-likelihood.md](../decisions/optimize-indel-contribution-to-likelihood.md)
