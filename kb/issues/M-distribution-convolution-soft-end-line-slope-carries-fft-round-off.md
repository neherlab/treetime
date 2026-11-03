# Soft-end straight line of a convolution takes its slope from FFT round-off

At an end of a convolution where an operand continues with a soft tail law, `fn reconstruct_neg_log_tails()` extends the neg-log values beyond the outermost resolved sample by a straight line ([kb/decisions/distribution-convolution-tails-by-exponential-tilting.md](../decisions/distribution-convolution-tails-by-exponential-tilting.md)). The slope of the line is `struct EndSlopes`: the difference of the outermost resolved sample and the sample 3 cells inward (`CONV_TAIL_MARGIN`), divided by 3 cells ([packages/treetime-distribution/src/distribution_ops/convolve.rs#L448-L463](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L448-L463), [#L589-L602](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L589-L602)). When the resolved range ends at the FFT trust floor ($10^{-13}$ of the peak, about 30 nats down), these samples carry the FFT round-off of about $10^{-3}$ nats. Divided by the width of 3 cells, the round-off becomes a slope error of order 0.1 nats per year, and the line carries it over its whole length. The round-off is not a smooth function of the operands: any change upstream, also one of $10^{-9}$ nats, draws a new round-off at the floor and thus a new slope.

When a posterior lies far out on such a line, its dates follow the slope error. The error does not shrink with the accuracy of the samples that the convolution resolves.

## Evidence

Nextstrain norovirus build of all genomes, v1 `timetree` with `--coalescent 0.01 --gen-per-year 50 --max-iter 2 --clock-filter 0 --clock-rate 0.00328 --confidence --reroot least-squares --branch-length-mode auto --divergence-units mutations-per-site --seed 1`, release build, one thread. Two runs of the same build differ only in how the exact tails at hard ends compute their gaps: by windowed tilts (the default), or by log-space quadrature for every gap. The windowed samples of the first run agree with exact quadrature over the full overlap within $3.6 \cdot 10^{-9}$ nats.

- **Backward pass**: the log marginal likelihood of the node times is $-16287.6179$ in both runs
- **First forward pass, node `NODE_0000228`** (a child of the root): the root posterior and the cavity (root posterior without the message of the child) agree within $1.1 \cdot 10^{-4}$ nats over the range of the posterior of the child. The forward message, the convolution of the cavity with the branch distribution, agrees within $2 \cdot 10^{-5}$ nats from 1965 to 2000. Its resolved range ends at 1962.24, 29.9 nats below its peak, and agrees there within $1.1 \cdot 10^{-3}$ nats. Older than that, the message is the straight line, with slopes of $-2.870$ and $-2.772$ nats per year in the two runs
- **Posterior**: the posterior of `NODE_0000228` lies from 1946.9 to 1955.0, 7 to 15 years beyond the start of the line and 50 to 74 nats down on it. The slope difference of 0.098 nats per year tilts the message by 0.79 nats across the posterior, and its mean moves from 1951.153 to 1951.198 years
- **Propagation**: the descendants of `NODE_0000228` inherit the shift in the forward pass, and the next iterations inherit it through the node times they use (polytomy resolution, the coalescent factor): the log marginal likelihood of the second backward pass differs by 0.074 nats, and the final dates differ by up to 0.42 years (`NODE_0000227`, 90% interval 4.8 years wide)

**Reproduction**: run the build above twice, once with a switch that makes `fn quadrature_is_cheaper()` of [packages/treetime-distribution/src/distribution_ops/exact_tails.rs](../../packages/treetime-distribution/src/distribution_ops/exact_tails.rs) return `true` for every gap, and compare the node dates by tip sets. To see the mechanism, write out `parent_time_dist`, `parent_except_subtree`, `dist_from_parent`, and `combined` of `NODE_0000228` in `fn refine_distribution_from_parent()` ([packages/treetime/src/timetree/inference/forward_pass.rs](../../packages/treetime/src/timetree/inference/forward_pass.rs)) in the first forward pass, and compare the slope of `dist_from_parent` older than 1962.

## Consequences

- **Dates on soft-end lines**: a node whose posterior lies on the straight line of a message gets dates that depend on the round-off at the trust floor of that message. On the build above, the dates move by up to 0.42 years between two runs whose samples agree within $10^{-8}$ nats
- **Accuracy of reference runs**: a reference run that computes every tail sample at hard ends by exact quadrature has the same lines at soft ends, so its dates carry the same error. On the build above, the reference and the run that computes every gap by quadrature differ by up to 0.10 years, and the default run by up to 0.32 years. A comparison of dates against such a reference cannot resolve differences below this size

> [!IMPORTANT]
> **Investigation required.** It is not shown why no tilt extends the older end of the forward message of `NODE_0000228`: its resolved range ends at the untilted trust floor, although the tilting loop of `fn neg_log_at_ends_by_tilting()` runs at a soft end too. Collect the slopes and the results of the tilts at that end, and the mass windows of the operands of `WindowedConvolution` that bound the convolution ([packages/treetime-distribution/src/distribution_ops/edge_convolution.rs](../../packages/treetime-distribution/src/distribution_ops/edge_convolution.rs)). Count how many posteriors of the build lie on a soft-end line, and how far out.

## Fix direction

> [!IMPORTANT]
> **Decision required.** The soft-end straight line is part of the approved tilting decision. Options, not yet chosen:
>
> - **Slope from accurate samples**: take the slope from samples whose round-off is at most $10^{-6}$ nats (for example above $10^{-7}$ of the peak of the FFT that produced them, or bounded by the round-off bound of `FftConvolver::overlap_sums_round_off()`), over a baseline of many cells, so that the round-off enters the slope with a small weight. The line stays an approximation of the tail
> - **Slope from the tail laws**: take the slope from the tail laws of the operands at that end, where they define one, instead of from sampled values
> - **Computed tails at soft ends**: extend the resolved range at a soft end inside the output window by the windowed tilts and the log-space quadrature of the exact tails, so that a posterior inside the window never lies on a line
>
> Evidence for the choice: the build above, before and after, compared between the two fill modes of the reproduction; the date spread between them measures the remaining sensitivity.

## Related

- [kb/decisions/distribution-convolution-tails-by-exponential-tilting.md](../decisions/distribution-convolution-tails-by-exponential-tilting.md): the tilts, the exact tails at hard ends, and the straight line at soft ends
- [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md): the tail laws of the messages
