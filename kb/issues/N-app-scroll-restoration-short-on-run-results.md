# Back navigation restores the run results page short of its scroll position

The router restores the window scroll position on back and forward navigation (`scrollRestoration: true` in `packages/app-ui/src/router.tsx`). On a run results page the restored position falls short: a page left at 1000 px came back at 767 px.

## Evidence

- TanStack Router restores the position in its `onRendered` handler, once the route has rendered
- The results page grows after that render: the tree workspace mounts behind `Suspense` in `packages/app-ui/src/runs/TreeView.tsx` and the Auspice tree lays out later, so at restore time the page is shorter than the saved position and the browser clamps the scroll offset
- Pages without asynchronous content restore exactly

## Required change

Reserve the final height of the tree workspace while it loads (a placeholder with the tree panel's height), so the page reaches its final height before the router restores the position.
