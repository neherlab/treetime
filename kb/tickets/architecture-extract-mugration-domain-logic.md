# Extract mugration domain logic from commands

## Description

`commands/mugration/run.rs` contains mugration business logic that belongs in a domain module. No `src/mugration/` module exists. The major domain extraction (GTR refinement, discrete marginal) was completed earlier, but mugration-specific computation remains in the command layer.

Domain functions in `commands/mugration/run.rs`:

- `execute_mugration` (86 lines): core mugration pipeline (state assembly, GTR construction, partition creation, marginal reconstruction, iterative refinement)
- `validate_weight_coverage` (20 lines): validates discrete attribute coverage against weights
- `compute_pi_from_weights` (10 lines): computes equilibrium frequencies from weight map
- `compute_pi_uniform` (3 lines): uniform equilibrium frequencies
- `apply_pseudo_counts` (9 lines): pseudo-count smoothing of equilibrium frequencies
- `WeightCoverageResult` struct

Domain types in `commands/mugration/output.rs`:

- `MugrationResult`: holds graph + partition + extracted traits + GTR output. Domain state container
- `MugrationTraitsOutput`, `MugrationGtrOutput`, `MugrationConfidenceOutput`: data extraction structs with `render_csv()` I/O methods

Create `src/mugration/` module:

- `src/mugration/mugration.rs`: `execute_mugration`, `MugrationInput`, `validate_weight_coverage`, `WeightCoverageResult`, `compute_pi_from_weights`, `compute_pi_uniform`, `apply_pseudo_counts`
- `src/mugration/result.rs`: `MugrationResult`, `MugrationTraitsOutput`, `MugrationGtrOutput`, `MugrationConfidenceOutput` (data extraction, not file I/O). The `render_csv()` methods can stay on the types since they produce strings rather than writing files.

Keep in `commands/mugration/`:

- `run.rs`: `run_mugration` (file I/O orchestration), `parse_mugration_input` (CLI arg parsing), `write_annotated_tree`, `write_gtr_json_file`, `write_confidence_csv`
- `comment_provider.rs` (Nexus output formatting)
- `input.rs` (CLI-specific input struct with arg-derived fields)

## Validation

- `src/mugration/` compiles independently with no `commands/` imports
- All existing tests pass (including `commands/mugration/__tests__/`)
- `commands/mugration/run.rs` imports from `crate::mugration` instead of defining functions locally

## Related issues

Source: [H-core-multi-client-architecture-library-purity](../issues/H-core-multi-client-architecture-library-purity.md)
