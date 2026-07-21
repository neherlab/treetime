# Mugration rejects tree leaves with no metadata entry

Mugration's `validate_trait_names()` ([`packages/treetime/src/partition/marginal_discrete.rs#L268-L274`](../../packages/treetime/src/partition/marginal_discrete.rs#L268)) hard-errors when a tree leaf has no matching row in the metadata file:

```
Mugration: tree leaves missing from metadata: <names>
```

This rejects valid datasets where some sampled tips lack a trait value, which is common when the trait column is sparsely populated.

## v0 behavior

v0 mugration ([`packages/legacy/treetime/treetime/wrappers.py#L777-L782`](../../packages/legacy/treetime/treetime/wrappers.py#L777)) iterates tree terminals and builds a pseudo-sequence per leaf, assigning `missing_char` to any leaf whose name is absent from the traits dict. Such leaves are treated as ambiguous (uniform profile) and reconstructed like any other missing observation; the run proceeds and prints the covered count ("Assigned discrete traits to N out of M taxa"). No error.

## Correct behavior

A tree leaf missing from the metadata should be assigned the missing-data character and the run should proceed, matching v0 and the FASTA subsystem, which fills absent tips with an ambiguous sequence ([`attach.rs#L61-L68`](../../packages/treetime/src/ancestral/attach.rs#L61)). The mugration model already carries a missing-data state (`DiscreteStates` is built with a missing symbol), so the assignment target exists.

Whether to also gate on very low coverage (as FASTA does at >1/3 missing) is a separate policy question, not required for v0 parity. The parity fix is: fill missing, do not error.

## Scope

The opposite direction (metadata names absent from the tree) was fixed to warn-and-proceed. This issue covers only the tree-leaf-missing direction, which still errors.

## Related

- [N-io-name-reconciliation-duplicated.md](N-io-name-reconciliation-duplicated.md) -- the same reconciliation is reimplemented per subsystem with inconsistent missing-data policy
- [M-io-sequence-name-matching-unreliable.md](M-io-sequence-name-matching-unreliable.md) -- broader name-matching reliability
- Design: [kb/proposals/input-name-matching-validation.md](../proposals/input-name-matching-validation.md) -- shared `reconcile_names` function; its cross-site table records this mugration missing-leaf hard-error as the strictest policy
