# App starts timetree runs that now stop for a missing sequence length

## Problem

`timetree` requires an alignment or `--sequence-length`, as v0 does (commit "feat(timetree): require an alignment or
--sequence-length, as v0 does"); without either, the run stops at the start with "--sequence-length is required
without an alignment". The app's run checks ([run_checks.rs](../../packages/app-backend/src/run_checks.rs)) only warn
that `--confidence` writes no intervals without an alignment, a sequence length or a clock rate std. dev.
(`confidence-without-rate-uncertainty`), so the app lets the user start a run that is bound to fail.

## Fix direction

Replace that warning by a check that blocks a timetree run without an alignment and without `sequence_length`, with
a fix that points to the sequence length setting. The `--confidence` part of the warning no longer applies: with a
sequence length the rate uncertainty is always available for an estimated rate.

## Validation

A timetree config without an alignment and without `sequence_length` is rejected before the run; with either, no
check fires.
