# Root Variable Indel Resolution: Default to Present

Variable indels at the root default to "present" (no gap). When the backward pass finds that some children have a gap and others have sequence at the same range, the root keeps sequence. The forward pass then resolves each child's indel state using the root (parent) state as tiebreaker.

This replaces the previous majority-rule behavior where the root adopted a gap when `deleted > present` among its children.

## Mechanism

The backward pass in `packages/treetime/src/seq/indel.rs` (`resolve_indels_backward`) classifies each range across children:

- All children gapped (or gapped + unknown + variable): resolved gap
- Mixed gap and sequence: variable indel (stored in `BTreeSet<(usize, usize)>`)
- No gap evidence: skipped

At the root, variable indels are not resolved. The `resolve_root_forward` function in `packages/treetime/src/ancestral/fitch_sub.rs` resolves substitution ambiguities but leaves indels untouched. The marginal-dense path in `packages/treetime/src/partition/marginal_dense.rs` clears `variable_indel` at the root without resolving.

For non-root nodes, `resolve_indels_forward` in `packages/treetime/src/seq/indel.rs` uses the parent state as tiebreaker: if the parent has a gap, the child inherits it; if the parent has sequence, the child keeps sequence. This is standard Fitch parsimony behavior and is unchanged.

## Previous behavior

The `Deletion` struct stored per-range counts (`deleted: usize`, `present: usize`) tracking how many children had gaps vs sequence. At the root, `resolve_root_forward` applied majority rule: if `deleted > present`, a gap was added. This required strict majority (ties defaulted to present).

## Trade-off

On binary roots (2 children), there is no behavioral change. A 1-vs-1 tie already defaulted to present under the old strict-majority rule.

On root polytomies (3+ children) where gaps outnumber sequence, the old behavior placed a gap at the root (fewer total indel events), while the new behavior keeps sequence (more deletion events on gapped children). Example: root with 3 children, 2 gapped - old behavior costs 1 insertion event, new behavior costs 2 deletion events.

Root polytomies with gap-majority variable indels are uncommon in practice. ML tree inference tools (IQ-TREE, RAxML, FastTree) produce fully resolved binary trees. Root polytomies appear in poorly resolved or consensus trees. Even when present, variable indels at the root are rare since indels are infrequent relative to substitutions in the viral datasets TreeTime targets.

## Rationale

Removing the per-child counts simplifies the data structure from `BTreeMap<(usize, usize), Deletion>` to `BTreeSet<(usize, usize)>`, reducing the indel tracking to presence/absence of disagreement. The forward pass already used parent state as tiebreaker for non-root nodes and never read the counts. The counts were only consumed at the root.

Defaulting to "present" at the root treats the ancestral sequence as complete, with gaps as derived states. This is a reasonable biological prior for rooted phylogenies where the root represents the most recent common ancestor.
