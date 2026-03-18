# Negative coalescent Tc accepted without validation

`--coalescent=-1.0` is accepted by clap and runs without error, producing
19/37 annotations (tips only, same pattern as
[internal node dates missing at scale](M-timetree-internal-dates-missing-scale.md)).
v0 crashes with NaN/log errors from the negative Tc.

Tc must be positive (it represents effective population size scaled by generation
time).
