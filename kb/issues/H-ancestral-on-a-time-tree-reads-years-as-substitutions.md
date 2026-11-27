# Ancestral on a time tree reads years as substitutions per site

`treetime ancestral` with one partition and an inferred model always sets the substitution rate to 1 after the fit, so the branch lengths of the input tree are read as expected substitutions per site. A time tree carries branch lengths in years, and no setting tells `ancestral` the scale of its input, so the reconstruction treats a branch of 85 years as 85 substitutions per site.

The rule is `InferredRate::for_partition_count()` ([packages/treetime/src/partition/inferred_rate.rs#L11-L13](../../packages/treetime/src/partition/inferred_rate.rs#L11-L13)): one partition gives `InferredRate::Unit`, whose `apply()` calls `gtr.set_mu(1.0)`. `fn build_alignment_partitions()` ([packages/treetime/src/ancestral/plan.rs#L180](../../packages/treetime/src/ancestral/plan.rs#L180)) applies it to every partition of the run. The rule follows v0 and is recorded under "Rates after the first fit" in [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md); the defect is that no input carries a different scale.

## The bundled pipelines run this case

`timetree` writes its Newick and Nexus trees with branch lengths in years, and has no divergence-unit tree output. The bundled Nextstrain pipelines pass that tree to `ancestral`:

```yaml
  - name: "ancestral"
    ancestral:
      tree: "{{ steps.timetree.outputs.nwk }}"
```

This appears in [data/dengue/2000/pipeline-nextstrain.yaml](../../data/dengue/2000/pipeline-nextstrain.yaml), [data/ebola/362/pipeline-nextstrain.yaml](../../data/ebola/362/pipeline-nextstrain.yaml), and the rsv, mpox and flu Nextstrain pipelines.

Measured on `data/ebola/362` with `treetime pipeline --config data/ebola/362/pipeline-nextstrain.yaml`, against `treetime ancestral` on the input tree of the same run with the same flags:

| Run | Tree total branch length | Reconstructed substitutions |
| --- | --- | --- |
| `ancestral` on the input divergence tree | 0.0542 substitutions per site | 1228 |
| `ancestral` on the time tree, as the pipeline runs it | 58.36 years | 3986 |

The pipeline reconstructs 3.2 times as many substitutions. On `data/dengue/20` the same comparison gives 4277 against 38621, a factor of 9. The inflation grows with the ratio of the time scale to the divergence scale, so it is largest for slowly evolving genomes and long time spans. Augur avoids the case: `augur refine` writes its tree in divergence units and `augur ancestral` reads that tree.

## No way to supply the rate of an inferred model

A named model keeps the `mu` of `--model-params`, and that path gives the right answer: `--model hky85 --model-params kappa=4 mu=0.000619` on the `data/dengue/20` time tree reconstructs 4257 substitutions, against 4277 on the divergence tree. The inferred model has no such path:

- `--model infer --model-params mu=...` is an argument error: "`--model-params` needs a named `--model`"
- No command reads a model document. `--output-model` writes `struct ModelOutput` ([packages/treetime/src/gtr/get_gtr.rs#L41](../../packages/treetime/src/gtr/get_gtr.rs#L41)), and only the app backend parses it, to display results. [kb/proposals/model-serialization-scalar.md](../proposals/model-serialization-scalar.md) proposes `--model-file` and is not approved
- The `timetree` model document reports `mu = 1.0` for an inferred partition. Its clock rate lives in a separate document (`timetree.clock-model.json`, `clock_rate`), so the substitution model and the time scale are never in one place

A pipeline that fits a model once and reuses it downstream, with `mu` in substitutions per year on a time tree, is therefore not expressible.

> [!IMPORTANT]
> **Decision required.** Three parts need a choice, and the first is enough to stop the wrong results.
>
> - **Guard**: `ancestral` could detect branch lengths that cannot be substitutions per site and stop with an error, warn, or stay silent. A tree whose total length far exceeds the number of sites, or whose root-to-tip spread matches the metadata dates, is evidence of a time tree. A threshold risks rejecting a legitimately saturated tree
> - **Scale input**: either let `--model-params mu=...` apply to `--model infer` (the model is fitted, the rate is given), add a `--branch-length-unit years|substitutions` that takes the rate from a clock model, or require the named-model path and document it
> - **Model reuse**: approve [kb/proposals/model-serialization-scalar.md](../proposals/model-serialization-scalar.md) so that `--output-model` of one command feeds `--model-file` of the next, and decide whether the document carries `mu` in substitutions per site, per year, or both with the unit named. [kb/issues/N-optimize-partition-rate-scale-invariance-tested-for-one-rate-step.md](N-optimize-partition-rate-scale-invariance-tested-for-one-rate-step.md) and the W normalization of [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md) ("A shared W is shared up to its scale") both bear on what the document must store

A fitted W is not reusable across partition schemes without care: the shared W of `--partition-by codon-positions` on `data/dengue/2000` has a C-T to A-G ratio of 3.12, and of `--partition-by rate-categories` 1.27, from the same pooled substitution counts, because the dwell-time term weights by the partition rates.

## Validation

- `ancestral` on a time tree and on the divergence tree of the same `timetree` run reconstruct the same substitutions, once the scale reaches the model
- The bundled Nextstrain pipelines reconstruct the substitution counts of an `ancestral` run on their input divergence tree
- A model written by `--output-model` and read back reproduces the run that wrote it, bit for bit, on the same tree and alignment

## Related

- [kb/proposals/model-serialization-scalar.md](../proposals/model-serialization-scalar.md) - the `--model-file` input that model reuse needs
- [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md) - the rule that sets mu to 1 after the first fit
- [kb/issues/N-optimize-reports-no-log-likelihood-or-trace-log.md](N-optimize-reports-no-log-likelihood-or-trace-log.md) - the missing likelihood that would show the cost of a wrong scale
