# Backward time messages are far too heavy far below their peak

A backward message of the time inference holds the exact convolution only near its peak. Below the trust floor of the FFT convolution, `CONV_TRUST_FRACTION` $= 10^{-13}$ of the peak [packages/treetime-distribution/src/distribution_ops/convolve.rs#L18](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L18), the message is a tail law reconstructed from a two-point secant (`fn reconstruct_neg_log_tails()` [packages/treetime-distribution/src/distribution_ops/convolve.rs#L355](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L355)) and then refitted from the outermost grid points. Where the data put a node's posterior far below the peak of a message it multiplies, the posterior reads this extrapolated tail instead of the convolution. On conflicted data the tail is hundreds of nats too heavy, and node dates in that region are wrong with and without the forward-posterior trim.

## Evidence

`data/flu/h3n2/200` at `--clock-rate=3.058e-3 --keep-root --time-marginal always --seed=1` (the rate that a default run estimates):

- **Posteriors far below their messages**: along `NODE_0000160 -> 161 -> 162 -> 163 -> 165 -> 167`, the posteriors of the 1972 cluster lie at the hard edge 1972.997 of their subtree distributions, 44 nats below the subtree's peak (1962.71) and 91 nats below the peak of the backward message from `NODE_0000162` (1948.53). The branch distribution of that edge peaks at 13.97 years and lies 392 nats above its peak at length 0
- **Tail too heavy**: against the direct integral of the branch density times the subtree density, the stored message is 335 nats too heavy at the parent time 1972.995 relative to 1970
- **Impossible mass**: in the untrimmed forward pass, `NODE_0000162` holds 23% of its mass below 1972.647, where its parent `NODE_0000161` holds $6.4 \cdot 10^{-6}$. In exact arithmetic a child's mass below a time is at most its parent's
- **Dates not converged**: the untrimmed pipeline at 4800 and 9600 grid points dates the cluster nodes 3 to 8 years apart; only 32 of 196 internal nodes converge to 1% of their 90% width between these two references. A cut of the parent at its $10^{-9}$ quantile removes the mass that the child puts there ([kb/decisions/timetree-forward-posterior-trimming.md](../decisions/timetree-forward-posterior-trimming.md#why-trimming-is-safe-and-where-it-is-not)). In a version of the forward pass with the trim and a spacing cap without the resolution check, dates of this region moved by 8 to 11 years against the run without both (Kolmogorov distance 1.0); most of that move came from the coarse forward messages of the cap
- **Ill-conditioned**: capping the forward spacing of a single node without a subtree distribution moves `NODE_0000138` by 0.35 of its 90% width; the uncapped 300-point run lies 0.55 to 1.4 widths away from the references on the cluster nodes
- **Against `rust`**: at this fixed rate, `NODE_0000061`, `062`, `064`, `065`, `083` and `084` differ from the code before the trim by 0.7 to 1.4 years, in a region where the 4800- and 9600-point references differ by 0.6 to 3.0 years. With default flags (estimated rate) no node of `data/flu/h3n2/200` moves by more than 0.02 of its 90% width

Reproduction:

```bash
./dev/docker/run just r treetime timetree --tree=data/flu/h3n2/200/tree.nwk \
  --metadata=data/flu/h3n2/200/metadata.tsv --alignment=data/flu/h3n2/200/aln.fasta.xz \
  --clock-rate=3.058e-3 --keep-root --time-marginal=always --seed=1 \
  --output-all=tmp/flu-deep-tail --output-confidence-tsv=tmp/flu-deep-tail/confidence.tsv
```

## Impact

- Node dates and intervals of nodes whose posterior lies tens of nats below a message peak, which happens where subtrees disagree on the timing (conflicted data, a fixed rate far from the one the data supports)
- The loss bound of the forward-posterior trim assumes exact messages, so the trim can move such nodes by more than their width
- The log marginal likelihood integrates the same messages ([H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md](H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md))

> [!IMPORTANT]
> **Investigation required.** How often a posterior lies more than 30 nats below the peak of a message it multiplies, on the example datasets with default flags and on the validation subset. The flu case above is at a fixed rate; the default run of the same dataset is not affected.

## Fix direction

Keep the dynamic range of a backward message where its product needs it: for example compute the deep tail of the convolution from the operands' tail laws in closed form instead of extrapolating the FFT result, or convolve in log space where the plain-space FFT loses digits. The operand tails of the FFT ([M-distribution-convolution-drops-operand-tails.md](M-distribution-convolution-drops-operand-tails.md)) are a related loss on the other end of the operands.

## Validation

- On the reproduction, every converged node of the 4800-point reference lies within 0.05 of its 90% width of the 300-point result, trimmed and untrimmed
- Against the direct integral of a parent message at the dates of its posterior's mass, the stored message is within 1 nat
