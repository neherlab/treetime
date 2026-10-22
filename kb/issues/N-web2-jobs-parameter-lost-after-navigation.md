# The jobs URL parameter of web2 is lost after the first navigation

web2 reads the URL parameter `jobs` (the thread count of the threaded module, `1` for the single-threaded module) once when the page loads. The router then replaces the URL without the query string, for example on the redirect from `/` to `/new`, so a reload on any later page starts with the default of 8 threads. The user sees no sign that the setting changed.

## Evidence

- `parseJobs(globalThis.location.search)` runs once at startup [packages/app-web2/src/main.tsx#L25](../../packages/app-web2/src/main.tsx#L25)
- The index route redirects to `/new` without the search string [packages/app-ui/src/router.tsx#L49](../../packages/app-ui/src/router.tsx#L49); in a browser, `location.search` is empty after this redirect

## Fix direction

Keep `jobs` across navigation, for example by carrying the search string through the router's redirects and links, or store the chosen value for the session of the tab and show it in the app.

## Validation

- After a load with `?jobs=1` and a navigation, a reload still runs on the single-threaded module
