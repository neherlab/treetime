# Newick branch lengths lose precision through f32 round-trip

The `bio` crate (v2.3.0) Newick reader parses branch weights as `f32` (`bio::io::newick`, `parse::<f32>()` in `bio-2.3.0/src/io/newick.rs`). The treetime-io Newick reader (`packages/treetime-io/src/nwk.rs`) widens to `f64` via `nwk_edge.weight as f64`, but the precision is already lost.

For example, `0.005` in a Newick string arrives as `0.004999999888241291` (the `f32` representation widened to `f64`). This affects all branch lengths read from Newick files. The relative error is bounded by f32 machine epsilon (~1.2e-7), so it does not affect scientific results at typical phylogenetic precision, but it prevents exact round-trip comparison of branch-length values in tests.

Observed in the optimize node data tests (`packages/treetime/src/commands/optimize/__tests__/test_augur_node_data.rs`), where assertions use `max_relative = 1e-6` to accommodate the precision loss.

## Possible fixes

- Replace the `bio` crate Newick reader with a custom parser that reads weights as `f64` directly. This would also fix the confidence/support parsing gap and enable Extended Newick support
- Upstream a PR to `rust-bio` changing the weight type from `f32` to `f64`
- Accept the precision loss and document the `max_relative` tolerance convention for Newick round-trip tests

## Related issues

All trace to the `bio` crate Newick parser (`bio::io::newick`) as root cause:

- [N-timetree-node-data-confidence-not-emitted.md](N-timetree-node-data-confidence-not-emitted.md) - the parser does not expose branch support values
- [N-timetree-unnamed-root-after-reroot.md](N-timetree-unnamed-root-after-reroot.md) - `assign_node_names` only runs during Newick parsing, not after topology changes (parser-adjacent)
- [kb/decisions/graph-based-phylogenetic-representation.md](../decisions/graph-based-phylogenetic-representation.md) - the `bio` parser does not support Extended Newick
- [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md) - the current `bio`-crate-based read path
