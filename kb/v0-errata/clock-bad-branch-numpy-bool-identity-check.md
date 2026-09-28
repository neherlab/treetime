# Clock regression tests the bad-branch flag by identity, which rejects numpy booleans

## v0 location

`ClockTree.setup_TreeRegression()` (`#ClockTree`, `#setup_TreeRegression`) [packages/legacy/treetime/treetime/clock_tree.py#L275](../../packages/legacy/treetime/treetime/clock_tree.py#L275)

## Erratum

v0 stores the flag `bad_branch` with two different types. `ClockTree._assign_dates()` sets it for an undated internal node with `np.all(...)`, which returns `np.True_` or `np.False_` ([packages/legacy/treetime/treetime/clock_tree.py#L152](../../packages/legacy/treetime/treetime/clock_tree.py#L152)). `TreeAnc._prepare_nodes()` sets the same flag with the builtin `all(...)`, which returns a Python `bool` ([packages/legacy/treetime/treetime/treeanc.py#L488](../../packages/legacy/treetime/treetime/treeanc.py#L488)).

The tip value of the root-to-tip regression tests the flag by identity:

```python
tip_value = lambda x: np.mean(x.raw_date_constraint) if (x.is_terminal() and (x.bad_branch is False)) else None
```

`np.False_ is False` is false, so a leaf whose flag is `np.False_` gets the tip value `None`, as if it were an outlier. The outlier mask of `TreeRegression.clock_plot()` tests the same flag by truthiness ([packages/legacy/treetime/treetime/treeregression.py#L509](../../packages/legacy/treetime/treetime/treeregression.py#L509)) and does not mask the leaf. The plot then calls `np.max` on an array that contains `None` ([L524](../../packages/legacy/treetime/treetime/treeregression.py#L524)) and fails with `TypeError: '>=' not supported between instances of 'float' and 'NoneType'`.

A leaf gets `np.False_` when rerooting turns an undated internal node into a leaf. Biopython's `Tree.root_with_outgroup()` does this to a root with one child. This is the crash reported in [neherlab/treetime#959](https://github.com/neherlab/treetime/issues/959) for `treetime` (timetree) and `treetime clock`. v1 shares the extra leaf itself: [kb/issues/M-reroot-single-child-root-becomes-extra-leaf.md](../issues/M-reroot-single-child-root-becomes-extra-leaf.md).

Two more identity tests in `ClockTree` fail the same way for internal nodes: `init_date_constraints()` tests `is True` ([packages/legacy/treetime/treetime/clock_tree.py#L389](../../packages/legacy/treetime/treetime/clock_tree.py#L389)), and `convert_dates()` tests `is False` ([L972](../../packages/legacy/treetime/treetime/clock_tree.py#L972)), so the warning for a node dated later than today is not shown for an internal node whose flag is `np.False_`.

## Evidence

- The adjacent outlier mask in `clock_plot()` tests the same flag by truthiness and disagrees with the tip value on the same node
- `TreeAnc._prepare_nodes()` computes the same flag with the same rule as a Python `bool`, so the result of the identity test depends on which function set the flag last, not on its value
- The regression itself skips `None` tip values ([packages/legacy/treetime/treetime/treeregression.py#L260](../../packages/legacy/treetime/treetime/treeregression.py#L260)), so the node silently drops out of the fit while the plot fails on it

## v0 impact

- `treetime` and `treetime clock` exit with the `TypeError` above when the reroot moves the root of a tree whose root has one child. `--keep-root` avoids it
- The same test excludes such a leaf from the root-to-tip regression although it is not an outlier. Leaves of the input tree get Python `bool` flags, so only internal nodes that become leaves are affected

## v1 status

v1 stores the outlier state as a Rust `bool` ([packages/treetime/src/clock/clock_state.rs#L227](../../packages/treetime/src/clock/clock_state.rs#L227)), so the same flag cannot have two representations. `treetime clock --reroot-tips=255211` on the tree from the issue exits with 0; the extra leaf appears in the output (see the linked issue) but does not crash the command.
