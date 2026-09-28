# The run log stops updating after switching between run tabs

Switching between the tabs of a running run (Results, Settings, Log) leaves the Log tab with the lines received before the switch. New lines never arrive, and after the run finishes the log stays incomplete until a page reload. For a time-tree run, the Log tab showed 73 lines up to 0.8 s after switching tabs while the run finished, and 225 lines after a reload.

## Evidence

- Each tab is its own route in `packages/app-ui/src/router.tsx`, so a tab switch unmounts one `RunPage` and mounts another. For a moment the run event query has no observer
- `fn runEventsQueryOptions()` in `packages/app-ui/src/api/events.ts` streams the events through `experimental_streamedQuery`, whose stream consumes the query's abort signal
- `removeObserver()` in `@tanstack/query-core` 5.100.11 (`src/query.ts`) cancels a running fetch that consumed its signal when the last observer leaves, so the stream closes. The cache keeps the lines received so far
- The query sets `staleTime: "static"`, so the observer of the new tab takes the partial cache as final and does not open the stream again

## Impact

- A run followed through its tabs shows a truncated log, although the run finished

## Required change

Keep the run event stream open across tab switches, for example by keeping one `RunPage` mounted for all tabs of a run, or restart the stream from the last received event when an observer subscribes to a partial stream. Cover the chosen behavior with a query-level test that removes and re-adds the observer during a stream.
