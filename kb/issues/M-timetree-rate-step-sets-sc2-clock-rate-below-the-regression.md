# Under a coalescent prior the rate step sets the sc2/4500 clock rate 30% below the regression

With a coalescent prior and no fixed `--clock-rate`, each refinement round moves the clock rate toward the maximum of the log marginal likelihood of the node times ($\log Z$, `log_lh_time`) under the frozen prior of the round ([kb/decisions/timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md), `fn step_marginal_clock_rate()` in [packages/treetime/src/timetree/optimization/marginal_clock_rate.rs](../../packages/treetime/src/timetree/optimization/marginal_clock_rate.rs)). On `data/sc2/4500` the first round targets $5.6 \cdot 10^{-4}$ substitutions per site per year, 0.61 times the rate of the root-to-tip regression, and the run ends at $6.6 \cdot 10^{-4}$. The Nextstrain ncov workflow fixes the SARS-CoV-2 rate at $8 \cdot 10^{-4}$ (`clock_rate: 0.0008` in `defaults/parameters.yaml` of `nextstrain/ncov` `8e16d94`).

## Evidence

The example config `data/sc2/4500/timetree.yaml` (`clock_filter: 3`, `coalescent: skyline`, `confidence: true`, `seed: 42`, two rounds):

| Step | Clock rate |
| --- | --- |
| Root-to-tip regression after the clock filter (500 tips flagged, $R^2 = 0.985$) | $9.20 \cdot 10^{-4}$ |
| Rate step of round 1: target | $5.61 \cdot 10^{-4}$ |
| Rate step of round 1: new rate (2/3 of the way in log rate) | $6.61 \cdot 10^{-4}$ |
| Rate step of round 2 (no maximum between the probes, the rate stays) | $6.61 \cdot 10^{-4}$ |

- `--coalescent opt` in place of the skyline: target $5.6 \cdot 10^{-4}$, final rate $6.6 \cdot 10^{-4}$, so the skyline does not cause the low target
- `--keep-polytomies`: the same rates, so polytomy resolution does not cause it
- The same targets ($5.61 \cdot 10^{-4}$ to $5.62 \cdot 10^{-4}$) on `rust` builds `602a051bb`, `d30093982` and `911cc83ea` with the same settings, so the target predates the normalization of `log_lh_time` and the latest changes to the time passes
- Without a coalescent prior each round keeps the regression: with the default flags (clock filter 4, no coalescent, `--seed 42`) the run ends at $8.54 \cdot 10^{-4}$

## Impact

The clock alone dates the tips that the clock filter flags and the internal nodes above them ([M-clock-filter-outlier-clade-dated-after-the-samples.md](M-clock-filter-outlier-clade-dated-after-the-samples.md)). At the low rate these nodes move to later dates: with the example config, 105 internal nodes lie after the latest input date 2026.36, by up to 1.64 years. With the regression rate fixed by `--clock-rate 9.2e-4`, the same run dates no internal node after the latest input date.

> [!IMPORTANT]
> **Investigation required.** Whether $\log Z$ of this tree under the frozen prior of the run peaks near $5.6 \cdot 10^{-4}$, or whether the three-probe parabola of the first round misplaces the maximum. Evidence to collect: the profile of `fn marginal_log_lh_at_rate()` over the rate with the frozen prior of round 1, from $4 \cdot 10^{-4}$ to $1.2 \cdot 10^{-3}$. If the profile peaks there, the next question is which part of the data pulls the rate down: the many branches without mutations, the 500 flagged tips that send no message, or the prior. [M-timetree-low-signal-rate-bias-claims-not-reproduced.md](M-timetree-low-signal-rate-bias-claims-not-reproduced.md) lists the open claims about a slow-rate bias of $\log Z$ on simulated data.

## Reproduction

```bash
./dev/docker/run just r treetime timetree --config data/sc2/4500/timetree.yaml \
  --output-all tmp/sc2-rate-step -v
```

The log lines `Clock rate stepped toward the maximum of the marginal likelihood of the node times` give the target and the new rate of each round.

## Validation

The profile of $\log Z$ over the rate is recorded here, and either the rate step reaches its maximum or the cause of a maximum far below the regression rate is identified.
