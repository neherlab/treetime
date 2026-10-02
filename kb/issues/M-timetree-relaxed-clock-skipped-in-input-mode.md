# Relaxed clock is skipped with input branch lengths

`fn relax_clock` in the refinement round returns early when the total sequence length is zero ([packages/treetime/src/timetree/round.rs#L250-L256](../../packages/treetime/src/timetree/round.rs#L250-L256)). With `--branch-length-mode=input` the branch model has no sequence partition, so `--relax` has no effect, with only an info-level log message. v0 runs the relaxed clock whenever it is requested ([packages/legacy/treetime/treetime/treetime.py#L311-L343](../../packages/legacy/treetime/treetime/treetime.py#L311-L343)).

## Impact

`--relax` combined with `--branch-length-mode=input` silently does nothing.

## Open question

v0 computes the relaxed-clock objective from the branch lengths and one-mutation scale; decide how v1 obtains the one-mutation scale without an alignment (`--sequence-length`), or reject the combination with an error.
