# Partitions cannot be configured with site ranges or per-partition models

A run builds one partition per `--alignment` file, and all partitions take the one `--model` choice of the command. No input assigns a different substitution model to each file, and no input splits one file into several partitions by site ranges, for example genes or codon positions in one alignment.

## Current state

- Each file becomes one partition that covers all its sites, in the order of the `--alignment` list (`fn split_sources()`, [packages/treetime/src/partition/result.rs#L49](../../packages/treetime/src/partition/result.rs#L49)). [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md) records the rule
- Optimize and timetree build their partitions with `fn build_partition_set()` [packages/treetime/src/partition/set.rs#L20](../../packages/treetime/src/partition/set.rs#L20), which passes the same representation and model to `fn build_marginal_partition()` ([packages/treetime/src/partition/create.rs#L23](../../packages/treetime/src/partition/create.rs#L23)) for every file [packages/treetime/src/optimize/pipeline.rs#L65](../../packages/treetime/src/optimize/pipeline.rs#L65) [packages/treetime/src/timetree/pipeline.rs#L327](../../packages/treetime/src/timetree/pipeline.rs#L327)
- Ancestral calls `fn build_marginal_partition()` once per file with the same model [packages/treetime/src/ancestral/plan.rs#L158](../../packages/treetime/src/ancestral/plan.rs#L158)
- Prune builds sparse JC69 partitions with `fn build_sparse_reconstruction()`, because it counts parsimony substitutions only [packages/treetime/src/prune/pipeline.rs#L58](../../packages/treetime/src/prune/pipeline.rs#L58)

## Impact

Segments or genes that need different models, for example a named model for one file and an inferred model for another, cannot share a run. Codon-position partitions and genes inside one alignment need a manual split into separate files first.

## Required boundary

A partition configuration layer consumes complete partition descriptions: the source file, the site range, and the model of each partition. It owns the assignment of sites and models and hands the descriptions to the construction functions above, so the construction policy (representation, GTR resolution, total length) stays in one place. Configuration sources, such as a config file ([kb/proposals/config-file-multi-partition.md](../proposals/config-file-multi-partition.md)), produce descriptions and do not duplicate construction. Scientific defaults and automatic partitioning remain separate decisions.

## Related issues

- [N-core-command-module-shared-ops-entanglement.md](N-core-command-module-shared-ops-entanglement.md)
- [N-io-multi-segment-genome-input.md](N-io-multi-segment-genome-input.md)
- [N-representation-infer-dense-stub.md](N-representation-infer-dense-stub.md)
