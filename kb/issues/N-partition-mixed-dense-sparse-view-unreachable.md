# Mixed dense and sparse partitions are unreachable

The optimize path carries two partition vectors and presents them as one combined view, dense entries first, then sparse ([`packages/treetime/src/optimize/run_loop.rs#L78`](../../packages/treetime/src/optimize/run_loop.rs#L78)). The likelihood, indel, and initial-guess passes sum contributions across that view, so the code reads as if one run can hold both representations.

No run does. Partition construction returns exactly one representation ([`packages/treetime/src/partition/create.rs#L17`](../../packages/treetime/src/partition/create.rs#L17)), and the optimize pipeline fills one vector and leaves the other empty ([`packages/treetime/src/optimize/pipeline.rs#L116`](../../packages/treetime/src/optimize/pipeline.rs#L116)). The `--dense` option selects the representation for the whole run. Multi-partition ancestral reconstruction builds one partition per gene and reconstructs each one separately, so it never mixes representations in a shared view either.

## Open question

Is a mixed dense and sparse run intended to become reachable?

- Keep the capability: define what a mixed run means for the quantities the combined view sums, and give it a construction path and a test. Until then the summing code is unverifiable.
- Drop the capability: one run holds one representation, the two vectors collapse to one, and the combined-view indirection disappears.

Settling this question defines whether the dense-first-then-sparse combined-view summation must support a real mixed run or can collapse to a single-representation vector.

## Validation

- Whichever way the question is decided, a test covers the decided behavior: a mixed run that produces a defined result, or a construction path that cannot produce a mixed run.
