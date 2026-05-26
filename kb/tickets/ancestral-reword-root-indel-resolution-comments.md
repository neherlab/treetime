# Reword root indel resolution doc comment and test comment

The doc comment at `packages/treetime/src/ancestral/fitch_sub.rs:109-112:` says "without per-child counts we cannot determine majority-rule direction." This frames the behavior as a limitation of a refactoring rather than a design choice. The test comment at `packages/treetime/src/ancestral/__tests__/test_fitch_sub.rs:237-238:` says "because we no longer track per-child counts" - describes what changed rather than what the algorithm does.

## Task

Reword the doc comment as a positive statement of the design choice, e.g.: "Variable indels at the root default to present (no gap). Direction is resolved in the forward pass via parent state."

Reword the test comment similarly, removing refactoring history.

## Related issues

Source: `.memory/15-pr-714/synthesis1.md` finding M2
