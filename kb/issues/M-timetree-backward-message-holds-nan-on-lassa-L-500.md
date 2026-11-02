# A time convolution gives NaN when an operand has two modes with a deep valley between them

The convolution of two finite time distributions produces NaN when one operand has two modes separated by a valley deeper than about 30 nats. The time inference then stops with "Cannot re-peak a distribution on [a, b]: its values contain NaN", in the backward pass (wrapped in "When sending the time message backward along edge E") or in the forward pass. The operands are finite and valid, so the NaN is a numerical defect of the convolution, not an input that the run must reject.

## Mechanism

`fn reconstruct_neg_log_tails()` [packages/treetime-distribution/src/distribution_ops/convolve.rs#L366-L379](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L366-L379) converts the FFT result to the neg-log form. It trusts every sample between the first and the last sample above the floor `peak * CONV_TRUST_FRACTION` [packages/treetime-distribution/src/distribution_ops/convolve.rs#L19](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L19), with `CONV_TRUST_FRACTION` $= 10^{-13}$, and takes $-\ln$ of each of them. Between two modes, the convolution falls below the floor, and there the FFT result is round-off of about $10^{-16}$ times the peak, with either sign. The logarithm of a negative sample is NaN. Any valley deeper than $\ln(10^{13}) \approx 29.9$ nats inside the trusted range can produce it.

## Evidence

In every traced case, the count of NaN samples equals the count of negative FFT samples in the trusted range:

| Run | Pass | NaN samples (all negative) | Most negative sample |
| --- | --- | ---: | ---: |
| Nextstrain norovirus GII.4 genome, edge 662 | backward | 107 | $-3.1 \cdot 10^{-17}$ |
| Nextstrain norovirus example data, all/VP1, node 450 | forward | 292 | $-3.0 \cdot 10^{-15}$ |
| `data/lassa/L/500`, edge 712 | backward | 107 | $-5.9 \cdot 10^{-16}$ |

The operands of the norovirus cases have two modes:

- **GII.4 genome**: the subtree distribution has a hard right end at its peak and a second mode 1.97 nats higher, with a barrier of 119 nats between them
- **Example data, all/VP1**: the message from the parent without the subtree has a barrier of 39.8 nats between its two modes

The branch likelihood, the coalescent factor, the products, the division, and the date constraints of these cases all give finite operands; the first NaN appears in the convolution, before the soft-tail fit and the resampling.

## Affected runs

- **Nextstrain norovirus builds**, all at the fixed clock rates of their workflows, all completed by v0: GII.4 genome (edge 662), all genome (edge 4093), the gene builds all/RdRp (edge 2858) and all/VP1 (edge 10382) in round 1, and the example data build all/VP1 (forward pass, node 450). The mechanism is traced on GII.4 genome and on the example data; the other three stop with the same error
- **`data/lassa/L/500`** with the input root kept or with the clock filter off:

  | Settings | Result |
  | --- | --- |
  | `--keep-root --clock-rate=0.00078` | finishes |
  | `--keep-root --clock-rate=0.00078 --confidence` | finishes |
  | `--keep-root --clock-rate=0.00078 --no-indels` | finishes |
  | `--keep-root --clock-rate=0.00078 --coalescent=opt` | NaN on edge 712 |
  | `--keep-root --clock-rate=0.00078 --clock-filter=0` | NaN on edge 752 |
  | `--reroot=midpoint` | finishes |
  | `--reroot=midpoint --clock-rate=0.00078` | finishes |
  | `--reroot=midpoint --clock-rate=0.00078 --coalescent=opt` | finishes |
  | `--reroot=midpoint --clock-rate=0.00078 --confidence` | finishes |
  | `--reroot=midpoint --clock-rate=0.00078 --clock-filter=0` | NaN on edge 826 |
  | `--reroot=midpoint --coalescent=opt --confidence --clock-filter=0` (estimated rate) | NaN on edge 826 |

  The failures are the same for every seed tried (1 to 5 and 42), and with `--keep-polytomies`. `--branch-length-mode auto` chooses `input` here (longest input branch 0.218 substitutions per site). The input root is the root that IQ-TREE wrote; with it and the fixed rate, the clock filter marks 56 outliers, many with apparent dates near the year 1100, and with the clock filter off these tips enter the node times. With the kept root, the coalescent prior starts from one segment over [1184.9, 2022.7] with the coalescent time scale $T_c = 2351$, and the failure comes in the first backward pass. The `timetree/lassa/L/500/keep-root-fixed-clock` row of `dev/smoke.toml` and `data/lassa/L/500/pipeline-nextstrain.yaml` (midpoint root, fixed rate, clock filter off) fail this way

The defect is not a regression: the code dates from `0e84f567c` (2026-08-13, "convolve via FFT with cropped tail reconstruction"), and the nightly of `a8f1633` stops on the norovirus example data in the same way.

## Reproduction

`test_edge_convolution_bimodal_operand_with_deep_valley_keeps_both_modes` in [packages/treetime-distribution/src/distribution_ops/__tests__/test_edge_convolution.rs](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_edge_convolution.rs) convolves an equal-weight mixture of $\mathcal{N}(20, 1)$ and $\mathcal{N}(80, 1)$ with $\mathcal{N}(0, 1)$ through `fn convolve_across_edge()`, once with the soft side left (backward) and once right (forward), where $\mathcal{N}(\mu, \sigma^2)$ is a normal density with mean $\mu$ and variance $\sigma^2$. The exact result is the equal-weight mixture of $\mathcal{N}(20, 2)$ and $\mathcal{N}(80, 2)$, so its neg-log density is equal at 20 and 80 and rises by $1/4$ from 20 to 21. Both cases stop with "Cannot re-peak a distribution on [8, 91.9999999999997]: its values contain NaN". The test is ignored until the fix:

```bash
./dev/docker/run just test-rs -E 'test(deep_valley)' --run-ignored only
```

On a dataset of the repository:

```bash
treetime timetree --tree=data/lassa/L/500/tree.nwk --metadata=data/lassa/L/500/metadata.tsv \
  --metadata-id-columns=accession --alignment=data/lassa/L/500/aln.fasta.xz \
  --keep-root --clock-rate=0.00078 --coalescent=opt --seed=42 --output-all=<dir>
```

## Fix direction

> [!IMPORTANT]
> **Decision required.** How `fn reconstruct_neg_log_tails()` treats the samples at or below the floor inside the trusted range. Options:
>
> - **Clamp to the floor**: replace each such sample by the floor before the logarithm. The result stays finite and keeps both modes; the valley holds the floor density instead of its true, much smaller value
> - **Fill the valley in log space**: interpolate the neg-log values across the untrusted samples from their trusted neighbors. The density keeps falling into the valley; the shape of the fill needs a rule
> - **Direct quadrature for deep valleys**: convolve such operands with `DirectConvolve` of [packages/treetime-ops/src/convolution.rs](../../packages/treetime-ops/src/convolution.rs), which has no FFT round-off, at a higher cost
>
> Evidence: in the traced cases the valley holds less than $e^{-30}$ of the peak density, so each option changes results only through that mass and through the tail fits next to the valley. The change needs a measurement of the dates and the log marginal likelihood on the affected runs.

## Validation

- Both cases of the deep-valley test pass, and the test drops its `#[ignore]`
- The norovirus builds above and the failing `lassa/L/500` settings complete; the smoke row and the pipeline config stop declaring the failure

## Related

- [N-distribution-mixed-nan-policy-undecided.md](N-distribution-mixed-nan-policy-undecided.md): this NaN arises inside the convolution from finite operands, so a NaN check at construction would only stop the run earlier
