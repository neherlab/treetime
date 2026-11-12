# `clock` fails on identical sampling dates without naming the cause

When every sample has the same sampling date, the root-to-tip regression has no spread on the date axis, so it cannot estimate a clock rate. `treetime clock` then fails with a message about a negative clock rate and a possible lack of temporal signal, and points to `--clock-rate`. The message does not say that the dates have no spread, which is the cause and which the user can check in the metadata.

## Reproduction

At commit `a71416d91`: the tree of `data/ebola/20` with a metadata table that gives every one of its 20 samples the date `2014.5`.

```sh
just r treetime clock --tree tmp/clock-same-dates/tree.nwk --metadata tmp/clock-same-dates/metadata.tsv -O tmp/clock-same-dates/out
```

```text
Error:
   0: The clock fit of the outlier filter requires a positive clock rate, also with --allow-negative-rate. Run with --clock-filter=0 to skip the outlier filter
   1: Clock rate is negative for all root positions. The data may lack temporal signal. Please specify --clock-rate explicitly.
```

With `--clock-filter 0` the run fails with `Clock model estimation failed` and the same cause line.

## Fix direction

Check the spread of the sampling dates before the regression, in the clock model that every command and the clock preview share, and fail with a message that names the cause: all dated samples have the same date (or fewer than two distinct dates), so the dates cannot give a clock rate; set `--clock-rate` or give samples with different dates. Test: a clock run on identical dates fails with that message.

## Validation

- The reproduction above fails with a message that names the missing spread of the sampling dates
