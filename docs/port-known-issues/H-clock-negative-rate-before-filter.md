# Clock command fails with negative rate before filtering outliers

The `clock` command fails with "Clock rate is negative for all root positions"
on datasets that v0 handles by filtering outliers first.

## Repro

```bash
./dev/docker/run ./dev/dev r treetime -- clock \
  --tree=data/dengue/100/tree.nwk --dates=data/dengue/100/metadata.tsv \
  --outdir=tmp/clock/dengue/100 --name-column=genbank_accession
```

## v0 vs v1 Behavior

**v0** (succeeds):

- Performs initial root finding
- Filters 8 outliers that violate the clock
- Re-runs regression with filtered data
- Gets positive rate: 6.44e-4 with R^2=0.93

**v1** (fails):

1. `estimate_clock_model_with_prefilter()` tries to get pre-filter clock model
2. Calls `get_clock_model()` -> `find_best_root()`
3. `find_best_root()` requires `has_positive_clock_rate()` for ALL candidate positions
4. ALL positions have negative rate (~-8.4e-4 with outliers included)
5. Fails before ever reaching the outlier filter step

## Root Cause

Location: `packages/treetime/src/commands/clock/find_best_root/find_best_root.rs:122-127`

```rust
if !has_positive_clock_rate(&best_res.clock_set) {
  return make_error!(
    "Clock rate is negative for all root positions. ..."
  );
}
```

v1 requires a valid (positive rate) clock model BEFORE filtering outliers, but
some datasets only have valid clock AFTER filtering. This is a chicken-and-egg
problem.

v0's `TreeRegression.find_best_root()` has a `force_positive=True` parameter
that can be set to `False` to allow initial root finding with negative rates.

## Error Message

```
Clock rate is negative for all root positions.
The data may lack temporal signal. Please specify --clock-rate explicitly.
```

## Affected Datasets

- dengue/100 (100 samples, 8 outliers)

## Potential Fixes

1. Allow negative rates during initial root finding (like v0's `force_positive=False`)
2. Apply outlier filtering BEFORE requiring positive rate
3. Use chi-squared minimization without rate sign check for initial pass
