# The skyline underestimates Tc on simulated data

> [!IMPORTANT]
> **Investigation required.** The cause is unknown, and it is unknown whether v0 has the same bias. Run v0 and v1 with the skyline on the same simulated trees (a few replicates per schedule at alignment length $L = 10^4$) and compare the coalescent time scale $T_c$ with the truth by depth. If v0 shows the same bias, whether to change the estimator is a scientific decision. If only v1 shows it, it is a v1 defect, and the differences below are the first suspects.

## Problem

The timetree suite of [validation/README.md](../../validation/README.md) (variant `keep-skyline-fixed`: keep root, true clock rate, `--coalescent-skyline`) estimates $T_c$ about 30% below the simulated value at $L = 10^4$. The error has a systematic dip near the root.

## Differences between v0 and v1

The skyline objective of both versions adds a smoothing penalty on $z_i = \ln T_{c,i}$, the log time scale of segment $i$, with stiffness $s$ (default 2 in both):

$$
P_{\text{v1}} = \frac{s}{2} \sum_i (z_{i+1} - z_i)^2, \qquad P_{\text{v0}} = s \sum_i (z_{i+1} - z_i)^2
$$

where:

- $P_{\text{v1}}$, $P_{\text{v0}}$ -- smoothing penalty of v1 and v0

The differences:

- Smoothing strength: v1 [packages/treetime/src/coalescent/skyline.rs#L360](../../packages/treetime/src/coalescent/skyline.rs#L360), v0 [packages/legacy/treetime/treetime/merger_models.py#L312](../../packages/legacy/treetime/treetime/merger_models.py#L312). At the same default, v1 smooths half as strongly. [kb/decisions/coalescent-skyline-convex-log-tc.md](../decisions/coalescent-skyline-convex-log-tc.md) does not record the factor
- Shape of $T_c$: v1 is piecewise constant on equal-width segments [packages/treetime/src/coalescent/skyline.rs#L60](../../packages/treetime/src/coalescent/skyline.rs#L60). v0 interpolates between evenly spaced points [packages/legacy/treetime/treetime/merger_models.py#L305](../../packages/legacy/treetime/treetime/merger_models.py#L305)
- When the skyline is fitted: [M-timetree-skyline-timing-v0-divergence.md](M-timetree-skyline-timing-v0-divergence.md)
- Node times it is fitted to: v1 marginal [kb/decisions/timetree-marginal-only-time-inference.md](../decisions/timetree-marginal-only-time-inference.md), v0 joint by default
- Lineage counts of the coalescent prior, fixed before the loop: [kb/decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md)

## Reproduction

From `validation/`, after `pixi run ./generate_all.sh`:

```bash
pixi run python run.py --suite timetree --set core --filter 'keep-skyline-fixed/.*/L10000_'
pixi run python evaluate.py --suite timetree
```

The `tc_log_bias` column of `results/timetree.tsv` holds the mean log error of $T_c$ per run.

## Validation

v1 and v0 skyline estimates on the validation trees have the same bias, or an approved decision records the difference.

## Related issues

- [M-coalescent-tc-band-undercovers-on-simulated-data.md](M-coalescent-tc-band-undercovers-on-simulated-data.md)
