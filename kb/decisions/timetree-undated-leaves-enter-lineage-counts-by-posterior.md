# Undated leaves enter the coalescent prior like dated leaves, through their posterior

A leaf without a date (or a clock-filter outlier) is a bad branch: it sends no message to its parent. Under a
coalescent prior its lineage is treated like that of any other leaf:

- the expected lineage count $\bar k(t)$ counts it with the survival function of its posterior time
- its posterior carries its leaf factor, $\exp(-\int^{x} \bar k_{-v}(t)/T_c\,dt)$, the survival of its branch against the
  other lineages, as the posterior of a leaf with a date range does

It still sends no message to its parent, so the marginal likelihood of the node times ($\log Z$, which the marginal
clock rate step maximizes) does not depend on a prior on its sampling date. The coalescent is conditional on the
sampling times and gives no such prior. With a flat prior, the message would be $\int L(\mu(x-t_p))\,dx =
\frac{1}{\mu}\int L\,db$: flat in the parent time but a factor $1/\mu$ per undated leaf in $Z$, which pulls the rate step
toward slow rates unless the prior matches the real sampling window. The leaf factor would make the message depend
weakly on the parent time through $\bar k$; that information is dropped.

An internal node whose children are all undated keeps the message from its parent as its posterior, without its own
factor. With weight $1 - n_\text{children} < 0$ that factor favors a later merger, and without child messages nothing
bounds it from the right.

A leaf with a partial date (a month or a year) is not a bad branch and already entered both ways.

**Type**: Scientific divergence from v0. v0 leaves bad branches out of the lineage counts while their parent's merger
term counts them ([packages/legacy/treetime/treetime/merger_models.py#L102-L105](../../packages/legacy/treetime/treetime/merger_models.py#L102-L105),
[packages/legacy/treetime/treetime/merger_models.py#L239-L250](../../packages/legacy/treetime/treetime/merger_models.py#L239-L250)),
and gives an undated leaf the message from its parent as its posterior, without a coalescent factor
([packages/legacy/treetime/treetime/clock_tree.py#L889-L891](../../packages/legacy/treetime/treetime/clock_tree.py#L889-L891)).
Before this decision v1 counted a bad branch as surviving at every time.

**v1 location**:

- `fn NodePosterior::summarize()` [packages/treetime/src/timetree/inference/result.rs#L111](../../packages/treetime/src/timetree/inference/result.rs#L111) stores the survival function of every posterior, bad branches included
- `fn LineageExpectation::from_posterior()` [packages/treetime/src/coalescent/mean_field.rs#L299](../../packages/treetime/src/coalescent/mean_field.rs#L299) reads it for every node; only a node without a posterior counts as surviving at all times
- The init seed: `fn GaussianClock::clock()` [packages/treetime/src/clock/gaussian.rs#L112](../../packages/treetime/src/clock/gaussian.rs#L112) reports the Gaussian marginal of every node without an included tip date. For an undated leaf it is the parent marginal carried down the branch, mean $t_p + \nu/\eta$ and variance $\sigma_p^2 + \nu/\eta^2$. `fn LineageExpectation::from_gaussian()` [packages/treetime/src/coalescent/mean_field.rs#L315](../../packages/treetime/src/coalescent/mean_field.rs#L315) uses it, untruncated. The Gaussian model itself is unchanged; it leaves excluded tips out
- The leaf factor: `fn refine_distribution_from_parent()` [packages/treetime/src/timetree/inference/forward_pass.rs](../../packages/treetime/src/timetree/inference/forward_pass.rs) multiplies the posterior of a leaf without subtree distribution by `fn apply_leaf_coalescent_factor()` [packages/treetime/src/timetree/inference/coalescent_factor.rs](../../packages/treetime/src/timetree/inference/coalescent_factor.rs) and normalizes it

## Why

Counting an undated leaf at every time adds a lineage between its true sampling date and the present that does not
exist. On the messy validation data (3% of tips undated, `validation/degrade.py`) at 100 tips per sample time there are
about 1.2 times as many undated tips as tips at the most recent sample time, so near the present $\bar k$ roughly doubles
and the pair rate $\bar k(\bar k-1)/2$ grows about five-fold. The fitted $T_c$ absorbs it, and the marginal rate step,
which maximizes $\log Z$ under that prior, pulls the clock rate down. The variants that estimate the rate by regression
(`none-est`) were unbiased on the same trees; gaps, ambiguous bases, runs of N and month-only dates had no measurable
effect. Giving the 134 undated tips of one set (`exp/rep00`, $L=10^4$, $\mu=10^{-3}$) their true dates removed the bias.

Leaving the leaf out of the counts, as v0 does, removes its lineage below its parent's merger, where it does exist,
while the parent's merger term still counts it. The posterior is what the time inference knows about the sampling date,
so it is the expectation the mean-field prior asks for. Without the leaf factor that posterior is the parent message
convolved with the branch, which at low signal is wide and late; the factor is the part of the coalescent that
constrains it.

## Validation

`validation/data_n100_messy` (4000 tips), `keep-opt` and `keep-skyline` variants, 10 replicates of 3 schedules per
cell, no failed runs. Median over runs; "clean" is `data_n100` on the same trees.

| quantity | $L$, $\mu$ | clean | messy, before | lineage counts only | this decision |
| --- | --- | --- | --- | --- | --- |
| rate error, `opt-est` | $10^4$, $10^{-3}$ | $-0.4\%$ | $-11.8\%$ | $-0.7\%$ | $+0.2\%$ |
| rate error, `opt-est` | $10^3$, $10^{-2}$ | $-0.7\%$ | $-12.0\%$ | $-1.0\%$ | $-0.2\%$ |
| rate error, `opt-est` | $10^3$, $3\cdot10^{-3}$ | $+5.0\%$ | $-11.6\%$ | $+2.9\%$ | $+5.9\%$ |
| $\ln(T_c/T_c^\text{true})$, `opt-fixed` | $10^4$, $10^{-3}$ | $0.013$ | $0.87$ | $0.042$ | $-0.014$ |
| $\ln(T_c/T_c^\text{true})$, `opt-fixed` | $10^3$, $10^{-3}$ | $0.150$ | $1.03$ | $0.436$ | $0.102$ |
| $T_c$, `opt-fixed`, `const` ($T_c^\text{true}=10$) | $10^3$, $3\cdot10^{-4}$ | $16.5$ | $41$ | $27$ | $17.0$ |
| date RMSE of undated tips (years), `opt-est` | $10^3$, $3\cdot10^{-4}$ | | $0.74$ | $0.60$ | $0.18$ |
| date RMSE of undated tips (years), `opt-est` | $10^3$, $3\cdot10^{-3}$ | | $0.18$ | $0.17$ | $0.091$ |
| date RMSE of undated tips (years), `opt-est` | $10^4$, $10^{-3}$ | | $0.108$ | $0.106$ | $0.077$ |
| internal-node date RMSE (years), `opt-est` | $10^4$, $10^{-3}$ | $0.038$ | $0.055$ | $0.040$ | $0.040$ |

The skyline variants agree. The validation runs use `--clock-filter 0` (`validation/run.py`), so they contain undated
tips but no clock-filter outliers.

Output on data where every leaf is dated and the clock filter flags no tip is byte-identical: the 13 `timetree` smoke
cases with a coalescent prior on `ebola/20`, `flu/h3n2/20` and `zika/20`. A clock-filter outlier is a bad branch too,
so a run with flagged outliers changes: on `ebola/362` with `--coalescent=skyline`, where the clock filter flags 18
tips, the clock rate moves from $8.04 \cdot 10^{-4}$ to $7.93 \cdot 10^{-4}$ ($-1.4\%$), the skyline changes shape
(segment boundaries move, $T_c$ rises up to 27-fold in the oldest segments and falls to 0.02 of before in the
youngest), and node dates move by up to 0.10 years. This effect is not validated against a known
truth ([kb/issues/M-timetree-undated-leaf-coalescent-counts-unvalidated-with-clock-filter-outliers.md](../issues/M-timetree-undated-leaf-coalescent-counts-unvalidated-with-clock-filter-outliers.md)).

Applying the internal factor
to internal nodes with only undated children, tried first, moved two such nodes of `exp/rep00`, $L=10^3$,
$\mu=3\cdot10^{-3}$ from $+0.04$ and $+0.33$ to about $+1.45$ years, and failed two trees whose right tails no longer
decayed.

## Open

[kb/issues/M-coalescent-lineage-grid-extends-past-the-latest-sample-with-undated-tips.md](../issues/M-coalescent-lineage-grid-extends-past-the-latest-sample-with-undated-tips.md):
at low signal the posteriors of undated tips reach decades past the latest sample, and so does the lineage grid.
