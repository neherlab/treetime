# Rebuild complete coalescent event state

Prevent topology refinement from running coalescent inference with missing node events.

## Required changes

- Rebuild internal time distributions and `likely_time` after every topology change.
- Make event collection fallible and reject missing required node times.
- Validate the expected event count and root/leaf/internal roles before likelihood evaluation.
- Preserve the complete previous inference state if rebuilding fails.

## Validation

- Polytomy refinement with and without topology changes.
- Missing-time error injection and whole-state atomicity assertions.
- Independent direct coalescent likelihood comparison.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/H-timetree-coalescent-events-incomplete-after-topology-change.md](../issues/H-timetree-coalescent-events-incomplete-after-topology-change.md)
