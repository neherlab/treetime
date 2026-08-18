# Combined hard-approach law: `mass` drops the linear slope term that `eval` keeps

A `HardApproachLaw` in the `Combined` regime carries both a power-law exponent $b > 0$ and a linear slope $s \neq 0$. Its two evaluators model different functions:

- `eval` returns the full two-term neg-log law $y_\text{edge} - b\ln(|t - t_\text{hard}| / D) + s\,(t - t_\text{edge})$, including the slope.
- `mass` returns $p_\text{edge}\,D / (b + 1)$, the pure power-law integral with the slope term **dropped**, where $D = |t_\text{edge} - t_\text{hard}|$ is the sub-grid gap width and $p_\text{edge} = e^{-y_\text{edge}}$.

So for a `Combined` law, `mass` is not the integral of $e^{-\texttt{eval}}$. The two agree only when $s = 0$ (the `Divergent` regime), which is the only regime a direct fit produces. `Combined` laws arise solely from composition: multiplying two `Distribution::Function` values whose facing tails are both `HardApproach` adds the shape terms of a divergent and a finite message, yielding $b > 0$ and $s \neq 0$ together.

The dropped-slope approximation is documented as "negligible over the sub-grid gap, mirroring the fit". That justification holds for a fitted law (where $s = 0$ by construction) but is unverified for a composed law, where $s$ is a sum of message slopes and $sD$ has no established bound.

## Affected code

- [`packages/treetime-grid/src/hard_approach_law.rs#L178-L186`](../../packages/treetime-grid/src/hard_approach_law.rs#L178-L186): `mass` collapses `Divergent` and `Combined` into one arm, discarding the slope.
- [`packages/treetime-grid/src/hard_approach_law.rs#L143-L162`](../../packages/treetime-grid/src/hard_approach_law.rs#L143-L162): `eval` keeps the slope term for `Combined`.
- [`packages/treetime-distribution/src/distribution_ops/multiply.rs#L359-L360`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L359-L360): the composition that produces `Combined` laws.
- [`packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L304-L312`](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L304-L312): `tail_mass` calls `HardApproachLaw::mass`, so a composed tail's mass is the slope-dropped value.

## Impact

The mass of a composed hard-approach tail is approximate. This feeds mass-domain sizing for distribution arithmetic, so the stored grid window for a product distribution can be sized from a mass value that omits the slope contribution of the sub-grid gap. The error scales with $e^{sD}$ over the gap; for small $sD$ it is below the numerical contract, but no invariant guarantees $sD$ is small after repeated composition.

## Exact form (pending decision)

The slope-inclusive integral over the gap has a closed form. With $u = |t - t_\text{hard}|$,

$$m = p_\text{edge}\, e^{sD} \int_0^{D} (u/D)^b\, e^{-su}\,du = \frac{p_\text{edge}\, e^{sD}}{D^{b}} \int_0^{D} u^{b}\, e^{-su}\,du .$$

For $s > 0$ this is a lower incomplete gamma, $\int_0^{D} u^{b} e^{-su}\,du = s^{-(b+1)}\,\gamma(b+1,\,sD)$, available as `statrs::function::gamma`. As $s \to 0$ it converges to $D^{b+1}/(b+1)$, recovering the current form; $s < 0$ needs the same analytic continuation with attention to overflow.

## Decision required

Two options, one scientific choice:

- **Keep the approximation** (current behavior): document that `Combined::mass` drops the slope by design and bound the acceptable $sD$. `eval` and `mass` remain inconsistent for `Combined` laws.
- **Make `mass` exact**: implement the incomplete-gamma form so `mass` equals $\int e^{-\texttt{eval}}$ for every regime. This changes numerical mass output for composed distributions and requires an independent oracle plus approval, per the project's numerical-correctness consent rule.

The current implementation preserves prior behavior (the approximation) and locks it with a test case; it does not adopt either option as a resolution.
