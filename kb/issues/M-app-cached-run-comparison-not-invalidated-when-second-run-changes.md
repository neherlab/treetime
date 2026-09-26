# A cached run comparison is not invalidated when the second run changes

The UI caches the answer of `GET /api/runs/{id}/compare/{other}` and relies on the app-wide event stream to invalidate it. When run `other` changes, for example when it finishes, the event lists only paths of run `other`, and none of them covers the comparison path, which lies under run `id`. The comparison page keeps showing the answer from before the change.

## Evidence

- `fn run_stale_paths()` in `packages/app-commands/src/runs/app_events.rs` lists `/api/runs` (exact), `/api/runs/{other}` (subtree) and `/api/clade-in-runs` (exact) for a change to run `other`
- `fn staleCoversKey()` in `packages/app-ui/src/api/keys.ts` matches a stale path against the leading path segments of a query key, so `/api/runs/{other}` covers `["api", "runs", other, ...]` and not `["api", "runs", id, "compare", other]`
- `Comparison` in `packages/app-ui/src/runs/ComparePage.tsx` reads the comparison with `staleTime: Infinity`, so the entry is read again only when it is invalidated

## Reproduction

1. Start a long time-tree run `b` and open the comparison of a finished run `a` with `b`; the comparison has no estimates, because `b` has no results yet
2. Wait until `b` ends with `ok`: the run records of `a` and `b` update, the comparison keeps showing no estimates until the page is reloaded

## Impact

- The comparison page shows the comparison from before the second run ended: no estimates and no ancestors
- Outside clients that cache answers by the stale paths of the event stream have the same gap

## Open question

The server must list a path that covers the comparison. Options:

- List `/api/runs/{id}/compare/{other}` for every pair that contains run `other`: the server does not know which pairs clients cache, so it would list one path per other run
- Add a scope that matches a path segment anywhere, for example "every path that contains run `other`", and give the comparison a key that such a scope covers
- Move the comparison to a path under both runs or to a query form, for example `/api/compare?runs=a,b`, and list `/api/compare` (exact) for every run change

## Related

- [kb/decisions/app-request-contract-full-command-config.md](../decisions/app-request-contract-full-command-config.md): event streams, stale paths and the UI data layer
