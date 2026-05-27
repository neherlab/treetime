# Implement per-position ambiguity mask for ancestral reconstruction

Compute a per-position mask identifying alignment columns where all tips are ambiguous. Used to filter spurious mutations and as a field in ancestral node data JSON output.

## Purpose

1. Mutation filtering: `collect_mutations()` skips mutations at masked positions. Without this, inferred bases at all-ambiguous columns produce spurious mutations (the inferred base is arbitrary when no tip provides signal).
2. Output field: `"mask"` string of `"0"`/`"1"` characters in ancestral node data JSON. `"1"` = masked. Not consumed by `augur export v2` but needed for VCF reconstruction and format completeness.

## Algorithm

From augur's `create_mask()` (`augur/ancestral.py:139-183`):

For FASTA input: for each alignment column, check if every tip has the ambiguous character (`'N'` for nuc, `'X'` for AA). Uses original tip sequences (`reconstructed=False`), not ML-inferred.

```
mask = [true; alignment_length]
for each tip:
  for each position:
    if tip_sequence[position] not in [ambiguous, gap]:
      mask[position] = false
```

A position is masked iff no tip provides a non-ambiguous, non-gap character at that position.

For VCF input: mask non-variable positions (positions not in the VCF). Union with the FASTA ambiguity check.

## Changes

1. Add `fn create_mask(partition, graph) -> Vec<bool>` that iterates leaf node sequences and checks per-position ambiguity against `alphabet.unknown()` and `alphabet.gap()`.
2. Use the mask in mutation collection: skip `edge_subs()` entries at masked positions when building `nodes.*.muts`.
3. Convert mask to `"0"`/`"1"` string for the `mask` field in ancestral node data JSON.
4. Wire into `run_ancestral_reconstruction()` output path.

## Locations

- New: mask computation function (near `packages/treetime/src/ancestral/` or in the node data JSON output module)
- `packages/treetime/src/partition/marginal_sparse.rs` - leaf sequence access via `SparseNodePartition.seq.sequence`
- `packages/treetime/src/partition/marginal_dense.rs` - leaf sequence access via profiles
- `packages/treetime/src/commands/ancestral/run.rs` - wire into output

## Related issues

- Source: mask subsection of [../reports/augur-node-data-json.md](../reports/augur-node-data-json.md)
