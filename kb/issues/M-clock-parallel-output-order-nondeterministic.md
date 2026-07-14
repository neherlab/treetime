# Parallel clock CSV row order is nondeterministic

Repeated runs of the same binary record differences in `rerooted.clock.csv`. `fn gather_clock_regression_results()` pushes rows into `ArrayQueue` during parallel traversal and collects queue order directly; CSV serialization preserves that scheduling-dependent order [packages/treetime/src/clock/rtt.rs#L28-L80](../../packages/treetime/src/clock/rtt.rs#L28-L80). The separately observed SVG difference is tracked in [N-clock-svg-output-nondeterministic-untraced.md](N-clock-svg-output-nondeterministic-untraced.md).

## Potential solutions

- O1. Produce keyed rows in parallel and sort by the approved topology/reference order before serialization.
- O2. Preallocate indexed result slots and fill each node’s deterministic slot. This avoids a sort when the ordering index already exists.

## Recommendation

Use the existing topology/reference-order index to construct deterministic slots. The row order must match the command's selected topology order and remain stable across worker counts.

## Related issues

- [N-reroot-missing-v0-golden-master.md](N-reroot-missing-v0-golden-master.md)
- [N-clock-svg-output-nondeterministic-untraced.md](N-clock-svg-output-nondeterministic-untraced.md)
