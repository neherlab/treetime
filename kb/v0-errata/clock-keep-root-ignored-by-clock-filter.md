# Clock `--keep-root` still reroots when the clock filter is on

## v0 location

`estimate_clock_model()` calls `TreeTime.clock_filter()` with `reroot=params.reroot or 'least-squares'` [packages/legacy/treetime/treetime/wrappers.py#L966-L970](../../packages/legacy/treetime/treetime/wrappers.py#L966-L970). It checks `params.keep_root` only after the filter, for its own reroot step [packages/legacy/treetime/treetime/wrappers.py#L978-L991](../../packages/legacy/treetime/treetime/wrappers.py#L978-L991).

## Erratum

`--keep-root` and `--reroot` are in one mutually exclusive group, and `--reroot` keeps its default `best` when `--keep-root` is given [packages/legacy/treetime/treetime/argument_parser.py#L170-L179](../../packages/legacy/treetime/treetime/argument_parser.py#L170-L179). `TreeTime.clock_filter()` reroots whenever its `reroot` argument is set: once before the outlier search and once after it [packages/legacy/treetime/treetime/treetime.py#L486-L489](../../packages/legacy/treetime/treetime/treetime.py#L486-L489) [L511-L513](../../packages/legacy/treetime/treetime/treetime.py#L511-L513). The clock filter is on by default (`--clock-filter 4.0`) [packages/legacy/treetime/treetime/argument_parser.py#L154-L162](../../packages/legacy/treetime/treetime/argument_parser.py#L154-L162). So `treetime clock --keep-root` reroots the tree unless `--clock-filter 0` switches the filter off.

## Evidence

- The help text of `--keep-root` says "don't reroot the tree" [packages/legacy/treetime/treetime/argument_parser.py#L172-L179](../../packages/legacy/treetime/treetime/argument_parser.py#L172-L179)
- Adjacent wrappers pass no root under `--keep-root`: `timetree` sets `root = None if params.keep_root else params.reroot` [packages/legacy/treetime/treetime/wrappers.py#L480](../../packages/legacy/treetime/treetime/wrappers.py#L480), and `arg` does the same [packages/legacy/treetime/treetime/wrappers.py#L337](../../packages/legacy/treetime/treetime/wrappers.py#L337). `TreeTime.run()` passes this root to its clock filter, so the filter does not reroot [packages/legacy/treetime/treetime/treetime.py#L246-L258](../../packages/legacy/treetime/treetime/treetime.py#L246-L258)
- `treetime clock --tree data/flu/h3n2/20/tree.nwk --dates data/flu/h3n2/20/metadata.tsv --sequence-length 1400 --keep-root` writes `.output.newick` with a root that differs from the input root. With `--clock-filter 0` added, the output root is the input root

## v0 impact

- `treetime clock --keep-root` reports the root-to-tip regression for a root that the user did not choose
- The rerooted tree is written as `.output.newick` or `pruned.newick`, the names that v0 uses for runs that keep the root
- The reroot leaves input support values on the wrong splits ([clock-newick-support-fused-into-node-names.md](clock-newick-support-fused-into-node-names.md))

## v1 status

v1 `clock --keep-root` keeps the input root while the clock filter is on, as [kb/algo/clock.md](../algo/clock.md) states. On `data/flu/h3n2/20`, the root of `clock.nwk` from `treetime clock --keep-root` is the input root.
