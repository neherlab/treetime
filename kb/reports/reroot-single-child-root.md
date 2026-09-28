# Crash and extra leaf when rerooting a tree whose root has one child

Reported in [neherlab/treetime#959](https://github.com/neherlab/treetime/issues/959) for 0.11.x. Reproduced on `master` at [`3aaffc5f`](https://github.com/neherlab/treetime/tree/3aaffc5f600f16d613d16b7228e6d317354e5993). All source links below point to that commit.

## Summary

- **Symptom**: `treetime` (timetree) and `treetime clock` exit with `TypeError: '>=' not supported between instances of 'float' and 'NoneType'` while plotting the root-to-tip regression
- **Cause**: the input root has exactly one child. Rerooting with Biopython turns this root into an extra leaf without a date. Two checks of the flag `bad_branch` disagree on this leaf, so the plot receives `None` as a date
- **Effect beyond the crash**: the extra leaf is written to the output trees (`rerooted.newick` has 101 leaves for 100 samples)
- **Fix**: remove an undated one-child root before rerooting, and make all `bad_branch` checks agree. The fix proposed in the issue is not needed (see [Options](#options))

## Reproduction

Data from the issue: [`reroot-single-child-root/tree.nwk`](reroot-single-child-root/tree.nwk) (the Newick string from the issue body) and [`reroot-single-child-root/leaf_dates.tsv`](reroot-single-child-root/leaf_dates.tsv) (the attachment). Commands run from the repository root with treetime installed from the checkout (`pip install -e .`):

```bash
d=kb/reports/reroot-single-child-root
out=$(mktemp -d)

treetime --tree $d/tree.nwk --dates $d/leaf_dates.tsv --sequence-length 29903 --outdir $out/timetree
treetime clock --tree $d/tree.nwk --dates $d/leaf_dates.tsv --sequence-length 29903 --outdir $out/clock
treetime --tree $d/tree.nwk --dates $d/leaf_dates.tsv --sequence-length 29903 --keep-root --outdir $out/keep-root
```

Derived inputs for the related defects:

```bash
# The same tree without the one-child root (root at the sample MRCA i373494)
sed -E 's/^\((.*)\)i373496:0;$/\1;/; s/:3\.34e-05;$/;/' $d/tree.nwk > $out/tree_no_stem.nwk

# An undated one-child node u1 on the branch above i359805
sed -E 's/\((336749:0\.0001671,95162:0\.0012706)\)i359805:0\.0001338/((\1)i359805:0.0001)u1:0.0000338/' \
  $out/tree_no_stem.nwk > $out/tree_inner_node.nwk

# A date for u1
{ cat $d/leaf_dates.tsv; printf 'u1\t2020-03-01\n'; } > $out/dates_inner_node.tsv

treetime --tree $out/tree_no_stem.nwk --dates $d/leaf_dates.tsv --sequence-length 29903 --outdir $out/no-stem
treetime --tree $out/tree_inner_node.nwk --dates $d/leaf_dates.tsv --sequence-length 29903 --outdir $out/inner-node
treetime --tree $out/tree_inner_node.nwk --dates $out/dates_inner_node.tsv --sequence-length 29903 --outdir $out/inner-node-dated
treetime --tree $out/tree_inner_node.nwk --dates $out/dates_inner_node.tsv --sequence-length 29903 --keep-polytomies --outdir $out/inner-node-dated-kp
```

Results on `master`:

| Run                   | Exit | Result                                                                                                                 |
| --------------------- | ---- | ---------------------------------------------------------------------------------------------------------------------- |
| `timetree`            | 1    | `TypeError` in `clock_plot`. The traceback appears above the buffered progress output                                  |
| `clock`               | 1    | Same `TypeError`. `rerooted.newick` is written before the plot and has 101 leaves, including `i373496`                 |
| `keep-root`           | 0    | Root `i373496` kept and dated 2020-01-22, before the sample MRCA `i373494` (2020-02-04). Rate `9.029e-04`              |
| `no-stem`             | 0    | 100 leaves, rate `8.452e-04`                                                                                           |
| `inner-node`          | 0    | Rate `8.457e-04`, where `no-stem` gives `8.452e-04`                                                                    |
| `inner-node-dated`    | 0    | `u1` missing from `dates.tsv`. `i359805` dated 2020-05-20, the same as in `inner-node`: the date of `u1` has no effect |
| `inner-node-dated-kp` | 0    | `u1` dated 2020-03-01 as given, `i359805` moves to 2020-04-26                                                          |

## Input tree

A one-child ("unary") root is valid Newick. It marks a stem: an ancestor above the most recent common ancestor (MRCA) of the samples.

- All branch lengths in this tree are whole mutation counts divided by the genome length 29903, rounded to four significant digits. The stem from `i373496` to `i373494` is one mutation (`3.34e-05`)
- Leaves have numeric IDs. Internal nodes are named `iNNNNNN` with numbers up to 373496

This pattern fits a subtree cut from a large mutation-annotated tree, with the ancestor above the sample MRCA kept as the root. The tool that produced the tree is not known.

## Root cause

1. **Reroot creates a leaf.** `TreeTime.reroot()` [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/treetime.py#L546)] reroots with Biopython's `Tree.root_with_outgroup()`, either directly [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/treetime.py#L624)] or through `TreeRegression.optimal_reroot()` [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/treeregression.py#L463)]. Biopython removes the path child from the old root. It deletes an old root with one remaining child, but attaches an old root with any other number of children as a node of the new tree [[src](https://github.com/biopython/biopython/blob/d7e4b8b19399668b09442a5b35765d9186b5f665/Bio/Phylo/BaseTree.py#L880-L899)]. A one-child root has no children left at this point, so it becomes a leaf. Biopython 1.88 (2026-08-06) and its current `master` behave this way
2. **The flag has a numpy type.** `ClockTree._assign_dates()` sets `bad_branch` of an undated internal node with `np.all(...)` [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/clock_tree.py#L152)]. The value is `np.False_` for the old root, because its subtree has dated leaves
3. **Two checks disagree.** The regression's tip value uses an identity test, `x.bad_branch is False` [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/clock_tree.py#L275)]. This is false for `np.False_`, so the new leaf gets the tip value `None`. The plot's outlier mask uses truthiness [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/treeregression.py#L509)], so the leaf is not masked
4. **The plot fails.** `clock_plot()` builds an object array that contains `None` [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/treeregression.py#L507)], and `np.max` raises [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/treeregression.py#L524)]. The regression itself skips `None` tip values [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/treeregression.py#L260)], so only the plot fails. Both commands reach the plot: `timetree` [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/wrappers.py#L564)] and `clock` [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/wrappers.py#L1072)]

## Related defects

- **Other identity checks on `bad_branch`**: `init_date_constraints()` tests `is True` [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/clock_tree.py#L389)] and `convert_dates()` tests `is False` [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/clock_tree.py#L972)]. For internal nodes with `np.False_`, the warning "node is later than today, but it is not marked as BAD" is never shown. `TreeAnc._prepare_nodes()` computes the same flag as a Python `bool` [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/treeanc.py#L488)], so the type depends on which function set it last
- **Dated one-child nodes lose their date**: `resolve_polytomies()` removes every non-root one-child node [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/treetime.py#L705-L709)], also one with a date constraint. TreeTime supports dates on internal nodes [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/clock_tree.py#L378-L387)]. The date is dropped without a message (runs `inner-node-dated` and `inner-node-dated-kp`)
- **An undated one-child node inside the tree changes the timetree result**: rate `8.457e-04` with the node and `8.452e-04` without it (runs `inner-node` and `no-stem`). The difference appears in both joint (`--time-marginal never`) and marginal (`--time-marginal always`) runs. `treetime clock` gives the same rate and root for both trees, so the difference comes from the timetree iterations. The cause is not traced
- **Direct use of `TreeRegression.optimal_reroot()`** [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/treeregression.py#L416)] on a tree with a one-child root also creates the extra leaf (101 leaves on this tree). `TreeRegression` knows only the tip values of leaves, so it cannot tell a dated root from an undated one. Direct use also requires a `bad_branch` attribute on every node, which the class does not set
- **Unused field in the residual clock filter**: `residual_filter()` sets `exact_date` only when `type(node) is float` [[src](https://github.com/neherlab/treetime/blob/3aaffc5f600f16d613d16b7228e6d317354e5993/treetime/clock_filter_methods.py#L23)]. A tree node is never a float, so the value is always `None`. The outlier table drops the column, so there is no output effect

## Options

### Fix proposed in the issue

The issue locates the failing expression correctly and suggests a local fix: convert the tip values to `float64` (so `None` becomes `NaN`) and use `np.nanmax`/`np.nanmin`. This stops the crash. It is not recommended as the fix:

- The extra undated leaf stays in the tree. It is written to `rerooted.newick`, `rtt.csv`, and the timetree outputs, and it is plotted as a point with no date
- The crash is the only visible sign that the tree changed during rerooting. Hiding it in the plot leaves the wrong tree in place

### Remove every one-child node when the tree is loaded

Not recommended. With `--keep-root`, the stem is part of the requested analysis: TreeTime keeps `i373496` as root and dates it (run `keep-root`). Removing it at load time drops this output. Dated one-child nodes inside the tree would also lose their date constraints.

### Move the date of a removed node to its child or to a leaf

Not valid. A date belongs to one point on a lineage. The child is a later point, separated by a branch whose duration is unknown, and estimating this duration is the task of TreeTime. The error equals the branch duration. At the clock rate of this dataset, one mutation corresponds to

$$
\frac{1}{\mu L} = \frac{1}{8.45 \times 10^{-4} \cdot 29903} \approx 0.04 \text{ years} \approx 14 \text{ days}
$$

where $\mu$ is the clock rate in substitutions per site per year and $L$ is the sequence length. A branch without mutations still has a positive duration, so the transfer is never exact. A leaf is a separate, later sample, so moving an ancestor's date to a leaf merges two different observations.

### Change a dated one-child node into a zero-length leaf

This is the usual representation of a sampled ancestor. In TreeTime it is only approximate: the branch to the new leaf gets a positive duration, so the node is only loosely tied to the date, the coalescent model counts an extra lineage, and the regression adds tip slack to the new leaf. Not needed, because TreeTime already supports dates on internal nodes.

### Recommended

1. **Remove an undated one-child root in `TreeTime.reroot()`, before the reroot.** Rerooting replaces the root by definition, and Biopython already removes an old root with two children. An undated one-child root carries no sample, date, or sequence. Removing it is exact:
   - Substitution probabilities along a branch satisfy $P(t_1)P(t_2) = P(t_1 + t_2)$, so joining two segments leaves the likelihood unchanged
   - Root-to-tip distances and their variances add along a path, so the regression is unchanged
   - A dangling branch to a node without data contributes a factor of 1 to the likelihood

   On this dataset the result equals removing the stem by hand: `dates.tsv` is identical to run `no-stem`

2. **Keep a dated one-child root.** After the reroot it becomes a dated leaf, which TreeTime handles like any other sample
3. **Leave the root unchanged with `--keep-root`**, where the stem is part of the requested tree
4. **Store `bad_branch` as a Python `bool` and test it by truthiness** in all places
5. **Check after the reroot** that every new leaf carries a date, and raise an error otherwise. This reports any other way of creating a leaf without data at its source, not later in plotting
6. **Keep dated one-child nodes in `resolve_polytomies()`**, so their date constraints stay in the inference

The inner-node rate difference, direct use of `TreeRegression.optimal_reroot()`, and the unused `exact_date` field are separate defects and are not part of this fix.
