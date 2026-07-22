# Timetree rejects a non-positive inferred clock rate

## Deviation

When the root-to-tip regression yields a non-positive clock rate, the `timetree` command errors instead of producing a time tree. v0 proceeds and writes output (a time tree, molecular-clock file, and ancestral reconstruction) using the negative rate.

The `clock` command does **not** share this behavior: it warns and continues (see [clock command](#clock-command-contrast) below).

## Rationale

Time-scaled inference maps divergence to time through `time = (divergence - intercept) / rate`. A non-positive rate makes this mapping meaningless: node times invert or diverge, and every downstream step (coalescent model, marginal time inference, confidence intervals) operates on nonsensical values. v0 computes and writes these values regardless; the output looks complete but is scientifically invalid.

Refusing early, with an actionable message, is preferable to emitting an invalid time tree. The error suggests specifying a known rate with `--clock-rate`, checking sampling dates, and verifying the alignment.

## Scope

The rejection fires wherever a non-positive rate reaches time inference:

- `--keep-root`: rerooting is disabled, so a negative-rate input root is used directly.
- `--allow-negative-rate`: see [below](#allow-negative-rate-in-timetree).
- Default rerooting already errors earlier ("Clock rate is negative for all root positions") when no root position yields a positive rate, matching v0's `Rerooting failed!`.

## `--allow-negative-rate` in timetree

The flag relaxes the **root search**: it lets the optimizer consider root positions with a negative slope (approaching midpoint rooting) rather than restricting the search to positive-slope roots. It does **not** disable the final positive-rate requirement. If the selected root still yields a non-positive rate, timetree errors as above, because time inference cannot proceed.

This is intentional: the flag widens rooting freedom, but a usable time tree still needs a positive rate. Users who want to inspect a negative-rate regression should use the `clock` command.

## Clock command contrast

The `clock` command reports the root-to-tip regression without performing time inference. A negative slope is a legitimate (if temporally uninformative) regression result, so `clock` warns and continues, matching v0. This makes `--keep-root` and `--allow-negative-rate` usable for diagnostics on data with weak temporal signal.

## v0 Reference

- `packages/legacy/treetime/treetime/wrappers.py:984`: `myTree.reroot(params.reroot, force_positive=not params.allow_negative_rate)` applies the flag to the reroot search only; v0 then proceeds with whatever rate results.
- `packages/legacy/treetime/treetime/treeregression.py:437-439`: `optimal_reroot` raises `Rerooting failed!` only when no acceptable root is found under `force_positive`, which v1 mirrors for the default path.

## v1 Implementation

- `ClockModel::from_regression` (`packages/treetime/src/clock/clock_model.rs`) errors on a non-positive rate. All timetree boundaries route through it (`packages/treetime/src/timetree/pipeline.rs`, `packages/treetime/src/timetree/optimization/reroot.rs`, `packages/treetime/src/timetree/refinement.rs`).
- `ClockModel::from_regression_allow_negative` warns and builds the model; only the clock command uses it, via `ClockRerootResult::into_clock_model_allow_negative` (`packages/treetime/src/clock/pipeline.rs`).
</content>
