# Timetree config files and pipeline steps disable the clock filter

`treetime timetree` runs with the outlier clock filter off when its arguments come from a config file (`--config`) or from a `treetime pipeline` step. The same run on the command line uses `--clock-filter=3.0`. The config file does not have to mention `clock_filter`: a one-line config that sets an unrelated key is enough.

## Reproduction

```bash
printf 'skyline_n_points: 4\n' > tmp/tt.yaml
treetime -v timetree --config tmp/tt.yaml --tree data/ebola/20/tree.nwk --metadata data/ebola/20/metadata.tsv ...
# argument dump: "clock_filter": 0.0

treetime -v timetree --tree data/ebola/20/tree.nwk --metadata data/ebola/20/metadata.tsv ...
# argument dump: "clock_filter": 3.0
```

## Mechanism

The field has two different defaults:

- The clap default is `3.0`: `clap(long, default_value = "3.0")` on `clock_filter` in [packages/app-cli/src/commands/timetree/args.rs#L394-L395](../../packages/app-cli/src/commands/timetree/args.rs#L394-L395)
- The struct derives `SmartDefault`, and the field has no `#[default]` attribute, so `TreetimeTimetreeArgsRaw::default().clock_filter` is `0.0`

The two argument sources that do not go through clap both start from `T::default()`:

- `fn overlay_config()` builds the merged value from `serde_json::to_value(T::default())`, applies the file, and then copies back only the CLI values whose source is `ValueSource::CommandLine` ([packages/app-cli/src/cli/config.rs#L38](../../packages/app-cli/src/cli/config.rs#L38)). The clap default `3.0` is not on the command line, so it is replaced by `0.0`
- Pipeline steps deserialize the step payload with `serde_json::from_value` and `#[serde(default)]` ([packages/app-cli/src/cli/pipeline/types.rs#L121](../../packages/app-cli/src/cli/pipeline/types.rs#L121))

The core runs the filter only when `params.clock_filter > 0.0` ([packages/treetime/src/timetree/pipeline.rs#L201](../../packages/treetime/src/timetree/pipeline.rs#L201)), so tips that violate the clock stay in the regression and in the time inference.

The generated schema `packages/schemas/input-config-timetree.schema.json` publishes `"default": 0.0` for `clock_filter`, while its description text says `Default=3.0`.

The sibling `clock` command declares both defaults (`#[default = 3.0]` in [packages/app-cli/src/commands/clock/args.rs#L195-L197](../../packages/app-cli/src/commands/clock/args.rs#L195-L197)) and is not affected. Every other field of `TreetimeTimetreeArgsRaw` and of the structs it flattens has matching defaults.

## Impact

- Wrong scientific results with no warning: clock outliers are not marked as bad branches, and they pull the clock rate, the root position, and the node dates
- The committed configs `data/flu/h3n2/200/timetree.yaml` and the timetree steps of `data/*/pipeline.yaml` run without the filter
- Results depend on how the arguments were supplied, not only on their values

## Fix direction

- Add `#[default = 3.0]` to `clock_filter`, and regenerate the schemas
- Prevent the class of bug: use one default source per field (for example `default_value_t = Self::default().<field>`, as `max_iter` and `skyline_n_points` already do), or add a test that compares the clap defaults with `Default` for every command
