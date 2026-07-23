# Domain payloads own file-format adapters

Core ancestral and timetree payloads directly implement Newick and Graphviz traits [`packages/treetime/src/payload/ancestral.rs#L20-L101`](../../packages/treetime/src/payload/ancestral.rs#L20-L101) [`packages/treetime/src/payload/timetree.rs#L116-L266`](../../packages/treetime/src/payload/timetree.rs#L116-L266).

The coupling includes presentation policy: `NodeTimetree::nwk_comments()` chooses the `date` key and two-decimal rendering beside scientific state [`packages/treetime/src/payload/timetree.rs#L129`](../../packages/treetime/src/payload/timetree.rs#L129). Payloads therefore change for domain, Newick, or Graphviz reasons.

## Dependency direction

The stable scientific payload layer imports format interfaces from the I/O layer. Ancestral and timetree payload families repeat this pattern for node labels, branch lengths, comments, and Graphviz attributes. Format evolution therefore modifies core types even when domain state is unchanged.

The issue is behavioral rather than a derive-only dependency: payload implementations choose keys, precision, omission rules, and presentation attributes.

## Design axes

- O1. Implement format adapters in boundary modules over small read-only payload access traits.
- O2. Convert payloads into explicit format projection values, then serialize those values.

Both keep scientific storage independent from format policy. O1 avoids intermediate ownership where streaming access is sufficient; O2 makes projection snapshots directly testable when formats require substantial transformation. Existing output formats and approved precision behavior must remain unchanged.

## Validation

- Format fixtures preserve Newick comments, branch lengths, labels, and Graphviz attributes exactly.
- Payload crates no longer import format-specific traits after migration.
- Domain tests construct payloads without I/O dependencies.
- Precision changes remain governed by the existing decision and are not bundled with adapter movement.

## Related decisions

- [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md)
- [kb/decisions/timetree-nwk-date-two-decimal-precision.md](../decisions/timetree-nwk-date-two-decimal-precision.md)
