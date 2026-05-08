# Extract mugration domain logic from commands/ to top-level module

## Description

Move domain algorithms from `commands/mugration/` to top-level `src/mugration/`. These import only from `gtr/` and `representation/partition/discrete` - no args coupling.

## What to move

- `gtr_refinement.rs` (381 lines) - iterative GTR inference from discrete trait data
- `discrete_marginal.rs` (170 lines) - discrete-trait marginal reconstruction (backward + forward)
- `comment_provider.rs` (30 lines) - partition comment formatting
- `representation/discrete_states.rs` (136 lines) - used exclusively by mugration

## What stays in commands/mugration/

- `args.rs` - CLI argument structs
- `run.rs` - command orchestration
- `input.rs` - input loading
- `output.rs` - output writing

## Related issues

- Source: [M-core-remaining-architectural-debt-after-extraction](../issues/M-core-remaining-architectural-debt-after-extraction.md)
