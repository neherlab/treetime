# Negative coalescent Tc accepted without validation

`--coalescent=-1.0` is accepted by clap and runs without error, producing
19/37 annotations (tips only, matching the
[internal dates missing for coalescent modes](M-timetree-internal-dates-missing-coalescent.md)
pattern). v0 crashes with NaN/log errors from the negative Tc.

Tc must be positive (it represents effective population size scaled by generation
time).
