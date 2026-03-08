# Internal node dates missing in nexus for coalescent, skyline, and input-BL modes

`timetree.nexus` contains date annotations only on tips when using fixed-Tc
coalescent (`--coalescent=<N>`), skyline (`--coalescent-skyline`), or input
branch lengths (`--branch-length-mode=input`). Default mode correctly annotates
all nodes (tips + internal).

## Scope

flu/h3n2/20 default: 37/37 annotations (19 tips + 18 internal).
flu/h3n2/20 with `--coalescent=10`: 19/37 annotations (tips only).

The issue is universal across all tested Tc values (0.1, 10, 1000), skyline
mode, and input-BL mode. `--coalescent-opt` and `--time-marginal` modes are not
affected.

## Root cause

`NodeTimetree.nwk_comments()` at
[`timetree.rs#L129`](../../packages/treetime/src/representation/payload/timetree.rs#L129)
writes the date annotation only when `self.time` is `Some`.

For the affected modes, `self.time` is `None` on internal nodes because:

- **Coalescent (fixed Tc) and skyline:** the second `run_timetree()` call at
  [`run.rs#L150`](../../packages/treetime/src/commands/timetree/run.rs#L150) with
  the coalescent prior overwrites time distributions from the first pass. The
  forward pass at
  [`forward_pass.rs#L47`](../../packages/treetime/src/commands/timetree/inference/forward_pass.rs#L47)
  calls `node.payload.set_time(time)` using `likely_time()` from the time
  distribution. With coalescent contributions, some distributions produce `None`
  from `likely_time()`, resetting the node's `time` to `None`.

- **Input-BL:** partitions are set to `vec![]` at
  [`run.rs#L65`](../../packages/treetime/src/commands/timetree/run.rs#L65), so
  marginal reconstruction does not run. Internal node time distributions are not
  fully computed, leaving `likely_time()` as `None`.

The tracelog shows `lh_pos=NaN` for coalescent runs, confirming positional
likelihoods are undefined for these internal nodes.

## Repro

```bash
./dev/docker/run ./dev/dev r treetime -- timetree --clock-filter=0 --coalescent=10.0 \
  --tree=data/flu/h3n2/20/tree.nwk \
  --dates=data/flu/h3n2/20/metadata.tsv \
  --outdir=tmp/repro-coal-dates data/flu/h3n2/20/aln.fasta.xz
grep -oP '\[&[^\]]*\]' tmp/repro-coal-dates/timetree.nexus | wc -l
# 19 (expected: 37)
```

## Related issues

- [Coalescent contributions use TBP coordinates, backward pass uses calendar time](H-timetree-coalescent-coordinate-mismatch.md)
  explains why coalescent contributions produce empty distributions
- [Coalescent CI excludes internal nodes](M-timetree-coalescent-ci-excludes-internal.md)
  compound effect when combined with `--time-marginal=only-final`
- [Internal node dates missing at scale](M-timetree-internal-dates-missing-scale.md)
  similar symptom from a different cause (numerical instability at scale)
