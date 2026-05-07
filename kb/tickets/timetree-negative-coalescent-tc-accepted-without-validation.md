# Negative coalescent Tc accepted without validation

`--coalescent=-1.0` is accepted by clap and runs without error, producing
tip-only annotations on flu/h3n2/20.
v0 crashes with NaN/log errors from the negative Tc.

The tip-only output is useful diagnostically. A negative Tc does not just represent an invalid parameter value. It also destabilizes the internal-node time inference path enough that only leaf dates remain materialized in the output.

Tc must be positive (it represents effective population size scaled by generation
time).

## Related issues

- Source: [N-timetree-negative-coalescent-tc.md](../issues/N-timetree-negative-coalescent-tc.md) -- delete after full resolution
