# Segment graphs with gamma rates or partition schemes are untested in optimize and timetree

The only test that combines segment graphs with gamma rate categories and partition schemes is the `ancestral` pipeline test [packages/treetime/src/ancestral/**tests**/test_pipeline_segment_graph.rs#L70-L82](../../packages/treetime/src/ancestral/__tests__/test_pipeline_segment_graph.rs#L70-L82) (dense and sparse marginal with `SiteModel::GammaRates`, odd and even partitions, also with Fitch). `optimize --segment-trees` and `timetree --segment-trees` have one test each, with uniform rates: [packages/app-commands/src/commands/optimize/**tests**/test_segment_trees.rs](../../packages/app-commands/src/commands/optimize/__tests__/test_segment_trees.rs) and [packages/app-commands/src/commands/timetree/**tests**/test_segment_trees.rs](../../packages/app-commands/src/commands/timetree/__tests__/test_segment_trees.rs). No test runs `--shared-params` with segment trees.

The `ancestral` test with gamma rates found a crash that review had not found: the sparse rate mixture read substitutions on edges that the segment does not carry, fixed in commit `63aa76246`. The branch-length fits of `optimize` and `timetree` run the same mixture and partition code on the graph, so they can hold defects of the same class.

## Validation

- `optimize`: with every edge shared by all segments, the fitted lengths equal those of the multi-partition run on the one tree, for gamma rates, partition schemes, and `--shared-params` (the oracle of the existing uniform-rate test in [packages/treetime/src/optimize/**tests**/test_pipeline_segment_graph.rs](../../packages/treetime/src/optimize/__tests__/test_pipeline_segment_graph.rs))
- `timetree`: the same graph with gamma rates and with partition schemes dates every node, and the dates equal those of the tree run
