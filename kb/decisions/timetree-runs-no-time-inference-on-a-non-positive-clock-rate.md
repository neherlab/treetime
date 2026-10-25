# Timetree runs no time inference on a non-positive clock rate

## Deviation

When the root-to-tip regression gives a clock rate at or below zero, `timetree` never runs the time inference on that rate. On most paths it takes the substitution rate of the Gaussian timetree instead; on the others it stops with an error. v0 proceeds with `--keep-root` and writes output (a time tree, molecular-clock file, and ancestral reconstruction) with the negative rate.

The `clock` command does **not** share this behavior: it warns and continues (see [clock command](#clock-command-contrast) below).

## Rationale

Time-scaled inference maps divergence to time through `time = (divergence - intercept) / rate`. A non-positive rate makes this mapping meaningless: node times invert or diverge, and every downstream step (coalescent model, marginal time inference, confidence intervals) operates on nonsensical values. v0 computes and writes these values regardless; the output looks complete but is scientifically invalid.

A non-positive slope on data with weak temporal signal does not mean that the data have no rate. The Gaussian timetree with the root time free has the log-likelihood $\kappa \ln \eta - A \eta^2/2 + B \eta - C/2$, maximal at $\eta = (B + \sqrt{B^2 + 4 A \kappa})/(2A)$, which is positive whenever $A > 0$, that is when the tip dates vary along the tree. Starting from this rate, the time inference and the marginal-likelihood step of the rate run as on any other input. On the validation suite the four runs that stopped on a non-positive slope (`const` rep03 and rep09, $L = 10^3$, `--keep-root`) complete with rate errors 0.46, 0.03, 0.03 and 0.76 and node RMSEs of 0.18 to 0.66 years, in line with their sets (`f6ca696c1`). The same fallback applies below the divergence floor of the start rate ([timetree-start-rate-falls-back-below-a-divergence-floor.md](timetree-start-rate-falls-back-below-a-divergence-floor.md)).

## Scope

With an estimated clock rate (a fixed `--clock-rate` has no regression rate), traced through `fn estimate_initial_clock()` in `packages/treetime/src/timetree/pipeline.rs`, `fn reroot()` in `packages/treetime/src/timetree/pre_loop.rs`, and `fn update_clock_model()` in `packages/treetime/src/timetree/round.rs`:

- **First clock fit**, every root method: without `--allow-negative-rate`, the `least-squares` root search keeps only roots with a positive slope, and stops the run when no root has one ("Clock rate is negative for all root positions", an input error, as v0 `Rerooting failed!`). On the chosen root, a slope at or below the divergence floor, and so every slope at or below zero, falls back to the Gaussian rate
- **Pre-loop refits** (before and after the clock filter), for a rerooting method: with `least-squares`, the refit keeps the regression rate and stops the run on a slope at or below zero ("Estimated clock rate is non-positive", from `ClockModel::from_regression`). Without `--allow-negative-rate` the root search keeps only roots with a positive slope, so the stop needs that flag. With `min-dev`, `midpoint`, `--outgroup` and `oldest`, a slope at or below the divergence floor falls back to the Gaussian rate
- **Clock model update of each refinement round**: a slope at or below zero falls back to the Gaussian rate
- **`--keep-root`**: the input root is used directly, so only the first fit and the rounds apply, and both fall back
- **No Gaussian rate either**: when a fallback finds no positive Gaussian rate (the tip dates do not vary along the tree, or the dated part of the tree has no substitutions), the run stops with an error that names the fit and both causes

## `--allow-negative-rate` in timetree

The flag relaxes the **root search** of `least-squares`: it lets the search consider root positions with a negative slope (approaching midpoint rooting) rather than restricting the search to positive-slope roots. It does **not** disable the positive-rate requirement of the pre-loop refits: if the selected root still yields a non-positive rate, `timetree` stops, because time inference cannot proceed. The other root methods do not search for a root by its slope, so the flag does not change them.

This is intentional: the flag widens rooting freedom, but a usable time tree still needs a positive rate. Users who want to inspect a negative-rate regression should use the `clock` command.

## Clock command contrast

The `clock` command reports the root-to-tip regression without performing time inference. A negative slope is a legitimate (if temporally uninformative) regression result, so `clock` warns and continues, matching v0. This makes `--keep-root` and `--allow-negative-rate` usable for diagnostics on data with weak temporal signal.

## v0 Reference

- `packages/legacy/treetime/treetime/wrappers.py:984`: `myTree.reroot(params.reroot, force_positive=not params.allow_negative_rate)` belongs to the v0 `clock` command (`fn estimate_clock_model()`), the only v0 command with `--allow-negative-rate` (`packages/legacy/treetime/treetime/argument_parser.py:447-458`); it applies the flag to the reroot search only, and v0 then proceeds with whatever rate results.
- v0 timetree has no `--allow-negative-rate`: `fn TreeTime.run()` reroots with `fn TreeTime.reroot()`, whose default is `force_positive=True` (`packages/legacy/treetime/treetime/treetime.py:542`, calls at `:260` and `:279`), so its reroots always require a positive rate. The v1 timetree flag has no v0 counterpart; see [kb/issues/M-timetree-allow-negative-rate-has-no-v0-counterpart.md](../issues/M-timetree-allow-negative-rate-has-no-v0-counterpart.md).
- `packages/legacy/treetime/treetime/treeregression.py:437-439`: `optimal_reroot` raises `Rerooting failed!` only when no acceptable root is found under `force_positive`, which v1 mirrors for the `least-squares` search.
- v0 has no fallback rate: `ClockTree.get_clock_model()` uses the regression slope as it is (`packages/legacy/treetime/treetime/clock_tree.py#L295-L311`).

## v1 Implementation

- `fn clock_fit_or_gaussian_rate()` (`packages/treetime/src/timetree/init.rs`) owns the fallback: `fn ClockFitResult::into_clock_fit_or()` (`packages/treetime/src/clock/clock_regression.rs`) keeps the regression rate above the threshold of the fit and otherwise takes the rate of `fn GaussianTree::rate()` (`packages/treetime/src/clock/gaussian.rs`).
- `ClockModel::from_regression` (`packages/treetime/src/clock/clock_model.rs`) errors on a non-positive rate; in `timetree` only the `least-squares` pre-loop refits reach it, through `fn ClockFitResult::into_clock_fit()`.
- `ClockModel::from_regression_allow_negative` warns and builds the model; only the clock command uses it, via `ClockFitResult::into_clock_model_allow_negative` (`packages/treetime/src/clock/pipeline.rs`).
