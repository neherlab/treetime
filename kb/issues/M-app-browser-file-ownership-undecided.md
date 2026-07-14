# Browser file ownership is undecided

The browser file picker records only a selected file's name and size. The bytes are not uploaded, while request construction treats the name as a server-local path. A browser-selected file can therefore refer to a path the server cannot access.

## Evidence

- `function FileSlot()` stores `{ name, size }` and discards the selected `File` object [packages/app-ui/src/components/FileSlot.tsx#L13-L43](../../packages/app-ui/src/components/FileSlot.tsx#L13-L43).
- `function buildDataPath()` forwards the stored name as a command path [packages/app-ui/src/components/RunButton.tsx#L9-L17](../../packages/app-ui/src/components/RunButton.tsx#L9-L17).

## Design axes

### Ownership boundary

- O1. Upload browser-selected bytes to a server-managed object with an opaque identifier. This supports arbitrary local files and makes lifetime, size limits, and cleanup explicit.
- O2. Remove browser file picking and permit only paths returned by the server dataset API. This keeps filesystem ownership entirely on the server but cannot accept arbitrary browser-local files.

### Lifetime

- O1. Scope uploaded data to a command and delete it when that command reaches a terminal state.
- O2. Scope uploaded data to a workspace with explicit deletion and expiry. This supports reuse but requires durable ownership and quota contracts.

## Recommendation

Use command-scoped uploads when arbitrary local input is a product requirement; otherwise expose only server-discovered datasets and remove the misleading browser picker. The product requirement must be approved before an implementation ticket can be created.

## Related issues

- [H-app-ui-displays-synthetic-results-and-ignores-inputs.md](H-app-ui-displays-synthetic-results-and-ignores-inputs.md)
