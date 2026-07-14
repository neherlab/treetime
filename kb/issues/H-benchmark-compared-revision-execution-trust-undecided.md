# Benchmark compared-revision execution trust is undecided

The benchmark harness extracts and runs the compared revision’s `dev/docker/run`. A revision under measurement can therefore define the mechanism intended to contain its own build execution.

`fn prepare_binary()` extracts the compared revision and invokes the extracted wrapper [dev/bench-graph-pass-cli#L151-L172](../../dev/bench-graph-pass-cli#L151-L172). This makes the measured revision part of the containment mechanism.

## Decision axes

### Containment owner

- O1. Use a separately identified trusted harness revision to build and execute both measured revisions.
- O2. Execute each revision’s harness. This measures historical tooling but grants the compared revision control of containment.

Recommendation: O1 for code benchmarking. Historical harness behavior can be inspected separately without execution.

## Recommendation

Use a separately identified trusted harness revision to build both measured revisions. Keep this issue ticketless until the containment-owner policy is approved.

## Related issues

- [H-benchmark-paths-cross-shell-filesystem-and-markdown-boundaries.md](H-benchmark-paths-cross-shell-filesystem-and-markdown-boundaries.md)
