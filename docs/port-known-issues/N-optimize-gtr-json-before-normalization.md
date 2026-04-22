# Optimize writes gtr.json before rate normalization

`write_gtr_json()` is called before `normalize_partition_rates()` in the optimize command. When the GTR model is inferred (`--model=infer`), normalization divides each partition's `gtr.mu` by the weighted-average rate and rescales all branch lengths. The JSON file records the pre-normalization `mu` value, so the serialized model does not match the in-memory model used for subsequent optimization or the output tree.

## Impact

The `gtr.json` output file contains a `mu` value that is inconsistent with the branch lengths in the output Newick/Nexus files. Downstream tools that read both outputs would use mismatched rate and tree parameters. The in-memory computation is unaffected because normalization mutates the partition's `gtr.mu` in place after the JSON write.

For named models (JC69, K80, etc.) where `mu` is set by the user, `normalize_partition_rates` is not called and the JSON is correct.

## Affected code

- Sparse write: [packages/treetime/src/commands/optimize/run.rs#L145](../../packages/treetime/src/commands/optimize/run.rs#L145)
- Dense write: [packages/treetime/src/commands/optimize/run.rs#L184](../../packages/treetime/src/commands/optimize/run.rs#L184)
- Normalization: [packages/treetime/src/commands/optimize/run.rs#L198](../../packages/treetime/src/commands/optimize/run.rs#L198)

## Fix

Move `write_gtr_json()` calls to after `normalize_partition_rates()`, or re-write the JSON after normalization.

## Related

- [M-optimize-gtr-mu-coordinate-mismatch](M-optimize-gtr-mu-coordinate-mismatch.md) - broader mu convention mismatch between optimizer and marginal propagation
