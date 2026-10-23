# Polytomy resolution fails with an internal error when a child of a polytomy is dated before its parent

`timetree` with polytomy resolution (the default) stops with an internal error in some runs, depending on the seed of the stochastic polytomy resolution:

```
Error:
   0: When running round 2
   1: Polytomy resolution failed
   2: Polytomy plan child 0 time must be finite and not older than its parent at 2.004783e3, got 2.004783e3. This is an internal error. Please report it to developers.
```

## Reproduction

```bash
./dev/docker/run just r treetime timetree --tree=data/flu/h3n2/500/tree.nwk --metadata=data/flu/h3n2/500/metadata.tsv \
  --alignment=data/flu/h3n2/500/aln.fasta.xz --seed=2251488387996782098 --output-all=tmp/timetree/flu/h3n2/500
```

- The x86_64 glibc dist binaries of `rust` `926cdf77b` fail with this seed at 1 and at 20 threads. The x86_64 musl dist binary passes with the same seed, so the outcome depends on the last bits of the math functions of the C library
- With random seeds, 1 of 400 runs of the glibc binary on this dataset failed; seeds 1 to 300 all passed
- The failing run first logs that 1 branch runs backwards in time and that its clock branch length was committed as zero

## Location

- `fn validate_plan()` in [packages/treetime/src/timetree/optimization/polytomy/apply.rs#L43-L60](../../packages/treetime/src/timetree/optimization/polytomy/apply.rs#L43-L60) rejects a child whose time is older than the time of the polytomy node (`child.time < parent_time`), without a tolerance
- The node times come from `state.time_inference.node_times()` in [packages/treetime/src/timetree/round.rs#L143-L162](../../packages/treetime/src/timetree/round.rs#L143-L162). A child with an observed or exact date keeps that date when it conflicts with the fitted clock, so it can be older than its parent

> [!IMPORTANT]
> **Investigation required.** Confirm that the rejected child is the branch that runs backwards in time (an observed leaf date or an exact date older than the parent), or a rounding difference between equal times (both print as `2.004783e3`). Then decide how polytomy resolution treats such a child: clamp its time to the parent, keep it out of the mergers, or reject the plan before validation with a user-facing message.
