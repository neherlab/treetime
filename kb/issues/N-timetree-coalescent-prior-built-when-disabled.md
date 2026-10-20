# The coalescent prior is built in runs without a coalescent prior

`fn run_refinement_loop()` builds a `CoalescentPrior` in every iteration and passes it on only when the coalescent mode wants a prior ([packages/treetime/src/timetree/refinement_loop.rs#L53-L59](../../packages/treetime/src/timetree/refinement_loop.rs#L53-L59)); the final step of the pipeline does the same ([packages/treetime/src/timetree/pipeline.rs#L399-L400](../../packages/treetime/src/timetree/pipeline.rs#L399-L400)). With `CoalescentMode::Disabled` (`fn CoalescentMode::wants_prior()` [packages/treetime/src/timetree/coalescent_timescale.rs#L185-L187](../../packages/treetime/src/timetree/coalescent_timescale.rs#L185-L187)) the prior is built and dropped.

`fn CoalescentPrior::new()` [packages/treetime/src/coalescent/mean_field.rs#L42](../../packages/treetime/src/coalescent/mean_field.rs#L42) evaluates $T_c$ on every cell of the lineage-expectation grid and returns an error when $T_c$ is not finite and positive. A run without a prior can therefore stop on an error of a prior it does not use, and spends the evaluation in every iteration.

> [!IMPORTANT]
> **Investigation required.** The run time of the unused construction on `data/sc2/4500` without a coalescent flag, and whether the timescale of a run without a prior can hold a $T_c$ that the constructor rejects.

## Fix direction

Build the prior only when the mode wants one, for example by returning `Option<CoalescentPrior>` from one constructor that takes the mode. [M-timetree-coalescent-modes-represent-invalid-states.md](M-timetree-coalescent-modes-represent-invalid-states.md) describes the related representation of the modes.

## Validation

A run without a coalescent flag constructs no `CoalescentPrior` (a test with a timescale that the constructor would reject completes).
