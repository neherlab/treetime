# Migrate command tests during commands/ dissolution

## Description

When `commands/` moves from the `treetime` library to consumer crates (`app-api`, `app-cli`), the tests in `commands/__tests__/` must be triaged and migrated:

- Tests that exercise CLI arg construction, file I/O, or end-to-end command execution through `run_*` functions: move to the consumer crate that owns the handler (likely `app-api` or `app-cli`)
- Tests that exercise domain functions (algorithm correctness, data transformations): keep in the library, move to the domain module's `__tests__/` directory
- Tests that use nonexistent file paths to verify fail-fast validation ordering: these test handler behavior, not domain logic. Move with the handler.

Unit tests are not suitable for testing file I/O operations. Integration tests or smoke tests are the right tool for verifying that the handler reads files, calls the pipeline, and writes output correctly.

## Specific tests to triage

Confidence tests in `commands/timetree/output/__tests__/`:

- `test_confidence_combine.rs`: tests `timetree::confidence::combine_confidence` (domain) -> move to `timetree/__tests__/`
- `test_confidence_extract.rs`: tests `timetree::confidence::extract_confidence_intervals` (domain) -> move to `timetree/__tests__/`
- `test_confidence_rate.rs`: tests `timetree::confidence::determine_rate_std` (domain) -> move to `timetree/__tests__/`

Mugration tests in `commands/mugration/__tests__/`:

- `test_discrete_marginal.rs`: tests domain marginal reconstruction -> move to `mugration/__tests__/`
- `test_gm_mugration.rs`: golden master test calling `execute_mugration` -> move to `mugration/__tests__/`
- `test_run.rs`: calls `run_mugration` with real files -> move to consumer crate
- `test_comment_output.rs`: tests `DiscreteCommentProvider` (domain) -> move to `mugration/__tests__/` or `partition/__tests__/`
- `test_augur_node_data.rs`: tests augur JSON formatting -> stays with augur builder location

Ancestral tests in `commands/ancestral/__tests__/`:

- `test_smoke_sample_from_profile.rs`: calls `run_ancestral_reconstruction` with real files + validation guard test with fake paths -> move to consumer crate
- `test_smoke_gtr_iterations.rs`: calls `run_ancestral_reconstruction` with real files -> move to consumer crate
- `test_augur_node_data.rs`: tests augur JSON formatting -> stays with augur builder

Optimize tests in `commands/optimize/__tests__/`:

- `test_augur_node_data.rs`: tests augur JSON formatting -> stays with augur builder

Timetree tests in `commands/timetree/__tests__/`:

- `test_pipeline.rs`: constructs `TreetimeTimetreeArgs` and calls `run_timetree_estimation` -> move to consumer crate
- `test_augur_node_data.rs`: tests augur JSON formatting -> stays with augur builder
- `test_auspice.rs`: tests auspice JSON output -> move to `timetree/output/__tests__/` (tests domain builder)

Stray tests outside commands/ that import `commands::shared::*`:

- `optimize/__tests__/test_args.rs`: tests clap parsing of opt-method args -> move to consumer crate
- `seq/__tests__/test_gap_fill.rs`: tests clap parsing of gap fill args -> move to consumer crate

## Related issues

- Source: [H-core-command-module-shared-ops-entanglement.md](../issues/H-core-command-module-shared-ops-entanglement.md)
- Source: [H-core-multi-client-architecture-library-purity.md](../issues/H-core-multi-client-architecture-library-purity.md)
- Prerequisite: `architecture-move-commands-to-cli-crate.md` (the dissolution itself)
