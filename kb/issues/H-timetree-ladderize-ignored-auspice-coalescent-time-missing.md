# Timetree: ladderize ignored in Auspice JSON; node times missing with coalescent

Two separate output bugs, both present when running `timetree`.

---

## 1. Ladderize / topology order not applied to Auspice JSON

### Symptom

`--ladderize=ascending`, `--ladderize=descending`, and `--topology-order` have no
effect on the Auspice v2 JSON output. The tree is written in the original traversal
order regardless of the flag.

### Root cause

[`run.rs`](../../packages/treetime/src/commands/timetree/run.rs) builds a reordered
graph and passes it to `write_tree_outputs`:

```rust
let plan = resolved.topology_order.plan(&output.graph)?;
let ordered = plan.ordered_graph(&output.graph)?;
// ...
write_tree_outputs(&ordered, &resolved.tree_outputs, &providers, Some(&auspice_ctx))?;
```

For Newick and Nexus output formats,
[`write_tree_outputs`](../../packages/treetime-io/src/graph.rs) uses the `ordered`
graph it receives. For `TreeWriteKind::Auspice` however, it ignores the `graph`
argument entirely and delegates to the `AuspiceWriter` trait object:

```rust
TreeWriteKind::Auspice => {
    writer.write_auspice(path)?;
}
```

The `TimetreeAuspiceCtx` passed as `auspice_writer` holds a reference to the
original `&output.graph`, not to `&ordered`. The Auspice JSON is therefore always
written in the original node insertion order.

### Fix direction

Either pass the ordered graph into `TimetreeAuspiceCtx` (replacing `output.graph`),
or change `AuspiceWriter::write_auspice` to accept a graph reference so
`write_tree_outputs` can supply the ordered one.

---

## 2. Node times missing from output when coalescent is specified

### Symptom

When `--coalescent`, `--coalescent-opt`, or `--coalescent-skyline` is specified,
the `numdate`, `date`, `clock_length`, and `branch_length` fields are absent (or
zero) for internal nodes in the augur node data JSON and Auspice JSON. Leaf times
from date constraints are unaffected.

### Repro

```bash
treetime timetree \
  --tree data/zika/20/tree.nwk --dates data/zika/20/metadata.tsv \
  --outdir tmp/out data/zika/20/aln.fasta.xz \
  --coalescent=10.0
# numdate is null for all internal nodes in the augur node data JSON
```

### Root cause

The backward pass initialises each internal node's accumulated distribution with the
coalescent contribution before processing children
([`backward_pass.rs:44-50`](../../packages/treetime/src/timetree/inference/backward_pass.rs)):

```rust
let mut result: Option<Distribution> = if node.is_leaf {
    None
} else {
    coalescent_contribs
        .and_then(|contributions| contributions.get(&node.key))
        .map(|contrib| (**contrib).to_plain())
};
```

`contrib` is `Arc<Distribution<NegLog>>` wrapping a `DistributionFormula`. Calling
`.to_plain()` on a `Formula<NegLog>` produces a `Formula<Plain>` whose closure
evaluates `exp(-neglog_value)` on every call
([`distribution.rs:405-413`](../../packages/treetime-distribution/src/distribution_core/distribution.rs)).

The neg-log values in the coalescent formula are large for almost all of the domain
(the coalescent contribution is small everywhere except near the true node time). In
float64, `exp(-x)` underflows to `0.0` when `x ≳ 745`. As a result, when
`multiply_formula_function` discretises this formula and multiplies with the child
message, the product array is all zeros. Calling `.normalize()` on an all-zero
`Function` returns `Distribution::Empty`. `likely_time()` on `Empty` returns `None`,
so `node.payload.set_time` is never called and `node.time` remains `None`.

The same underflow would affect small-Tc scenarios because the integral merger rate
`I(t)` grows proportionally to `1/Tc`, making neg-log values even larger.

### Fix direction

The coalescent contribution should not be applied via `to_plain()` on the neg-log
formula. The refactored architecture proposed in
[`kb/algo/coalescent-contribution-refactor.md`](../algo/coalescent-contribution-refactor.md)
addresses this directly: after all child messages are accumulated into `result` (a
`Function<Plain>` on a known grid), the coalescent weight is applied in-place as a
pointwise multiply using a scalar-only `CoalescentModel::eval(t, multiplicity)`
method. This avoids the neg-log → plain conversion on the formula closure entirely
and never produces an all-zero intermediate.

A targeted short-term fix would be to apply the coalescent contribution in log-space
by accumulating all distributions as `NegLog` throughout the backward pass and
converting to `Plain` only at the end, avoiding the underflow at the formula
evaluation stage.
