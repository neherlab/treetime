# Warn instead of error when mugration metadata has names not in tree

`validate_trait_names()` in [`packages/treetime/src/partition/marginal_discrete.rs#L275-L281`](../../packages/treetime/src/partition/marginal_discrete.rs#L275) hard-errors when metadata contains names not present as tree leaves. This rejects valid datasets where the metadata file has more samples than the tree (common when samples are pruned upstream).

## Cross-command consistency contract

v1 already handles this scenario correctly in other subsystems. Mugration should follow the same pattern:

| Subsystem | Extra names (not in tree) | Model to follow                                                                                     |
| --------- | ------------------------- | --------------------------------------------------------------------------------------------------- |
| FASTA     | Silently ignored          | [`attach.rs#L21-L22`](../../packages/treetime/src/ancestral/attach.rs#L21)                          |
| Dates     | Warning, continues        | [`date_constraints.rs#L92-L113`](../../packages/treetime/src/clock/date_constraints.rs#L92)         |
| Mugration | **Hard error** (fix this) | [`marginal_discrete.rs#L275-L281`](../../packages/treetime/src/partition/marginal_discrete.rs#L275) |

The dates subsystem's `warn_unused_date_constraints()` is the best model: it collects unused names, logs a warning with the count and a sample (first 10), and continues.

## What to change

In `validate_trait_names()` (`marginal_discrete.rs#L275-L281`), replace the `make_error!` with a `warn!` call for the `missing_in_tree` check. Keep the `missing_in_metadata` check (tree leaves missing from metadata, lines 267-273) as a separate question -- v0 assigns `missing_char` to those leaves rather than erroring, but that is out of scope for this ticket.

Before:

```rust
let missing_in_tree: IndexSet<String> = trait_names.difference(&leaf_names).cloned().collect();
if !missing_in_tree.is_empty() {
  return make_error!(
    "Mugration: metadata names missing from tree leaves: {}",
    missing_in_tree.iter().join(", ")
  );
}
```

After:

```rust
let missing_in_tree: IndexSet<String> = trait_names.difference(&leaf_names).cloned().collect();
if !missing_in_tree.is_empty() {
  let sample = missing_in_tree.iter().take(10).join(", ");
  let suffix = if missing_in_tree.len() > 10 { "..." } else { "" };
  warn!(
    "Mugration: {} metadata names not present in tree: {sample}{suffix}",
    missing_in_tree.len()
  );
}
```

## Test update

Update the existing test `test_discrete_marginal_attach_traits_rejects_metadata_name_missing_from_tree` in [`packages/treetime/src/commands/mugration/__tests__/test_discrete_marginal.rs#L55-L67`](../../packages/treetime/src/commands/mugration/__tests__/test_discrete_marginal.rs#L55). It currently asserts that extra metadata names produce an error. Change it to assert the function succeeds (returns `Ok`) and verify the warning is emitted (or just verify `Ok` and rename the test to reflect the new behavior).

## Verification

After the fix, the tb/20 mugration smoke test should pass:

```bash
./dev/docker/run ./dev/dev r treetime -- mugration --tree=data/tb/20/tree.nwk --states=data/tb/20/metadata.tsv --attribute=cluster --name-column=strain --outdir=tmp/verify/mugration/tb/20 --pc=1.0
```

Expected: succeeds with a warning about G22694 and G22721 not being in the tree.

## Related issues

- Source: [kb/issues/M-mugration-extra-metadata-names-rejected.md](../issues/M-mugration-extra-metadata-names-rejected.md) -- delete after full resolution
