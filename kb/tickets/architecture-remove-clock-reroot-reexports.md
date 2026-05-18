# Remove unused pub use re-exports from clock/reroot.rs

## Description

`clock/reroot.rs` re-exports `EdgeMergeInfo`, `EdgeSplitInfo`, `RerootChanges`, `RerootResult` from `partition/algo/topology_cleanup/reroot`. No production caller imports these through the re-export path. Creates dual import paths for `RerootChanges`.

## Fix

Remove the `pub use` line from `clock/reroot.rs`. Callers that need these types import directly from `partition/algo/topology_cleanup/reroot`.

1 line deletion.

## Related issues

- Source: [M-core-remaining-architectural-debt-after-extraction](../issues/M-core-remaining-architectural-debt-after-extraction.md)
