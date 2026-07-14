# Connect UI inputs, dataset selection, results, and errors

Make the application execute the visible configuration and render only the resulting typed scientific data.

## Required changes

- Move every command parameter into typed central state and build requests from it.
- Replace the complete file selection atomically when datasets change.
- Accept only server-discovered dataset paths; keep browser file selection unavailable until its ownership contract is approved.
- Persist typed command results and render them instead of mocks.
- Preserve and display typed errors with causes.
- Disable or remove synthetic result panels from completed real runs.

## Validation

- Parameter propagation for every command family.
- Dataset replacement and stale-path cases.
- Real success result, backend error, cancellation, and stale-response cases.
- End-to-end web and desktop parity.
- Full frontend and Rust checks.

## Related issues

- Source: [kb/issues/H-app-ui-displays-synthetic-results-and-ignores-inputs.md](../issues/H-app-ui-displays-synthetic-results-and-ignores-inputs.md)
