# Remove dead parameters from `resolve_root_forward`

`fn resolve_root_forward` at `packages/treetime/src/ancestral/fitch_sub.rs:113-125:` accepts five parameters but only reads three (`sequence`, `variable`, `chosen_state`). Both `gaps: &mut [(usize, usize)]` and `_variable_indel: &BTreeSet<(usize, usize)>` are never accessed in the function body. The `_variable_indel` prefix makes the disuse explicit, but `gaps` is silently dead.

The caller at `packages/treetime/src/ancestral/fitch.rs:218-224:` passes both. Tests at `packages/treetime/src/ancestral/__tests__/test_fitch_sub.rs:225:` and `:245:` construct data for both dead parameters.

## Task

Remove both dead parameters from the function signature. Update call sites in `fitch.rs` and `test_fitch_sub.rs`.

## Related issues

Source: `.memory/15-pr-714/synthesis1.md` finding M1
