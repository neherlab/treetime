# Node times one ulp apart count as a backward branch and stop polytomy resolution

An inferred node can be dated one ulp (unit in the last place) after the exact date of one of its children, although its posterior ends at that date. Two consumers of the node times compare them without the coordinate tolerance of the project and treat the pair as a branch that runs backwards in time:

- **Polytomy resolution** stops the run with an internal error: "Polytomy plan child 4 time must be finite and not older than its parent at 2.026325e3, got 2.026325e3. This is an internal error. Please report it to developers." The message prints both times with 6 digits, so they look equal
- **The time inference** warns after each pass that "N branch(es) run backwards in time" and commits their clock lengths as zero, a warning that names a conflict between an observed date and the clock

## Evidence

- **Polytomy error**: in a Nextstrain seasonal-flu build, the rejected child keeps its exact sampling date, and its parent lies exactly one ulp later, $2.27 \cdot 10^{-13}$ years at the year 2026. The run first warns about one backward branch, and that branch is the same edge
- **Warnings**: in four Nextstrain production builds, every branch that the warning counts is exactly one ulp backwards: 368 counts over the passes of lassa GPC, 38 of norovirus p22, 7 of measles N450, and 19 of the H5N1 cattle-outbreak HA build. The warning appeared in 13 of 21 replayed Nextstrain builds
- **Rounding dependence**: `data/flu/h3n2/500` with `--seed=2251488387996782098` stops with the polytomy error with the x86_64 glibc dist binaries of `rust` `926cdf77b` at 1 and 20 threads, and passes with the musl binary; 1 of 400 random seeds failed with glibc, and seeds 1 to 300 passed. Whether the one-ulp difference appears depends on the last bits of the math functions of the C library
- **Property test**: [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md) records a simulated case whose tip, dated exactly at 2024.3333333333333, lies $2 \cdot 10^{-13}$ years before its parent

```bash
./dev/docker/run just r treetime timetree --tree=data/flu/h3n2/500/tree.nwk --metadata=data/flu/h3n2/500/metadata.tsv \
  --alignment=data/flu/h3n2/500/aln.fasta.xz --seed=2251488387996782098 --output-all=tmp/timetree/flu/h3n2/500
```

## Mechanism

- **Node times**: `fn committed_time()` [packages/treetime/src/timetree/inference/forward_pass.rs#L291](../../packages/treetime/src/timetree/inference/forward_pass.rs#L291) commits `max(likely_time, parent_time)` for every inferred node, and a node with an exact date keeps that date ([kb/decisions/timetree-uncertain-leaf-dates-are-inferred.md](../decisions/timetree-uncertain-leaf-dates-are-inferred.md)). The posterior of a parent ends at the exact date of its youngest child, and its most likely time lies at that end when the branch carries no substitution

> [!IMPORTANT]
> **Investigation required.** The operation that places the end of the parent posterior, and with it the most likely time, one ulp after the exact date of the child: the grid arithmetic of the hard end, a resampling, or the conversion between dates and times before present.

- **Polytomy plan check**: `fn validate_plan()` rejects a child with `child.time < parent_time` [packages/treetime/src/timetree/optimization/polytomy/apply.rs#L52](../../packages/treetime/src/timetree/optimization/polytomy/apply.rs#L52)
- **Backward-branch count**: the commit of the clock lengths counts every edge with `duration < 0.0` [packages/treetime/src/timetree/inference/runner.rs#L210](../../packages/treetime/src/timetree/inference/runner.rs#L210) and warns from that count [packages/treetime/src/timetree/inference/runner.rs#L221-L234](../../packages/treetime/src/timetree/inference/runner.rs#L221-L234)
- **Owner of the comparison**: `fn coordinates_coincide()` [packages/treetime-grid/src/grid.rs#L188](../../packages/treetime-grid/src/grid.rs#L188) treats two time coordinates within 4 ulps of their scale as one point. The division of distributions already uses it for hard ends that rounding moves by one ulp ([M-timetree-outputs-depend-on-the-rounding-of-the-platform.md](M-timetree-outputs-depend-on-the-rounding-of-the-platform.md), section "Amplifying steps")

## Fix direction

Compare node times through `fn coordinates_coincide()` in both consumers, and in any other comparison of node times that tests their order, so that times one ulp apart count as equal. This is independent of the projection of inverted posterior modes, which [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md) leaves to a decision.

## Reproduction

`test_apply_plan_accepts_a_child_dated_one_ulp_before_its_parent` in [packages/treetime/src/timetree/optimization/polytomy/__tests__/test_apply.rs](../../packages/treetime/src/timetree/optimization/polytomy/__tests__/test_apply.rs) applies a plan without mergers to a polytomy whose parent is dated at 2026.3246575342469 and whose first child at the adjacent float 2026.3246575342466, the two times of the seasonal-flu build. It fails today with the internal error and is ignored until the fix:

```bash
./dev/docker/run just test-rs -E 'test(one_ulp_before_its_parent)' --run-ignored only
```

## Validation

- The test passes and drops its `#[ignore]`
- The flu command above finishes with the glibc build
- The four builds above no longer warn about backward branches that are one ulp long
