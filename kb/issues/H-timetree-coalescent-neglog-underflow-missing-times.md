# Timetree coalescent likelihood underflow removes internal node times

## Symptom

When `--coalescent`, `--coalescent-opt`, or `--coalescent-skyline` is specified, the `numdate`, `date`, `clock_length`, and `branch_length` fields can be absent or zero for internal nodes in the Augur node-data JSON and Auspice JSON. Leaf times from date constraints are unaffected.

## Reproduction

```bash
treetime timetree \
  --tree data/zika/20/tree.nwk --dates data/zika/20/metadata.tsv \
  --outdir tmp/out data/zika/20/aln.fasta.xz \
  --coalescent=10.0
# numdate is null for all internal nodes in the Augur node-data JSON
```

## Root cause

The backward pass initializes each internal node's accumulated distribution with the coalescent contribution before processing children in [`packages/treetime/src/timetree/inference/backward_pass.rs`](../../packages/treetime/src/timetree/inference/backward_pass.rs):

```rust
let mut result: Option<Distribution> = if node.is_leaf {
  None
} else {
  coalescent_contribs
    .and_then(|contributions| contributions.get(&node.key))
    .map(|contrib| (**contrib).to_plain())
};
```

The contribution is a `Distribution<NegLog>` wrapping a `DistributionFormula`. Calling `to_plain()` evaluates `exp(-neglog_value)` directly. When every neg-log value on the relevant grid exceeds the floating-point underflow threshold, all relative distinctions collapse to zero. Multiplication with the child message then produces an all-zero function, and normalization returns `Distribution::Empty`. `likely_time()` on an empty distribution returns `None`, so the internal node time is never assigned.

Small coalescence-time scales are particularly susceptible because the integral merger rate grows in proportion to the inverse coalescence time.

## Fix direction

Combine the coalescent contribution and the first child message in neg-log space, preserving the approved [coalescent-first multiplication order](../decisions/coalescent-multiplication-ordering.md). Before converting the combined distribution to plain probabilities, subtract its minimum neg-log value. This makes the peak probability one while preserving every likelihood ratio and prevents an all-zero intermediate.
