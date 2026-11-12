# Web2 keeps runs in the tab and exports them as run packages

**Type**: Application architecture without a v0 counterpart.

**Status**: Accepted and implemented.

## Decision

- **Runs live in the tab**: the runs of web2, their inputs, and their outputs live in the in-memory file system of the back-end worker of the tab ([app-web2-wasm-worker-back-end.md](app-web2-wasm-worker-back-end.md)). Nothing stores a run in browser storage. A reload, a closed tab, or a crash of the back-end worker ends every run of the tab, and each tab has its own runs
- **Browser storage holds the preferences and the form only**: `localStorage` keeps the preferences under `treetime-preferences` (the sidebar width, shared by every tab) and the theme under `theme`; `sessionStorage` keeps the unfinished form of each tab under `treetime-draft` (`packages/app-ui/src/preferences/storage.ts`). After a reload, the form keeps its settings and fetches its example inputs again; uploaded files belonged to the lost back end, and the user chooses them again
- **The run package keeps and shares a run**: a run leaves the tab as a run package, one zip with the record, the log, the inputs, the outputs, the YAML configuration, and the command line (`packages/app-backend/src/runs/package.rs`). Importing the package in any tab or in the desktop app recreates the run with equal results, and the CLI reproduces it from the extracted folder
- **Desktop keeps its runs folder**: the desktop app stores each run in its runs folder with the same back-end code (`RunStore`), so its history outlives the app. The two histories differ in lifetime: desktop lists every run of the runs folder, web2 lists the runs of one tab
- **The interface states the model**: web2 asks before the tab unloads while it holds runs, shows a banner when the back-end worker crashes, explains a run link without a run in the tab, offers "Import a run package" there, shows the size of each run, and deletes runs to free tab memory

## Options weighed

- **Runs that persist in the browser** (rejected): the back-end worker mirrors the run store into the origin private file system (OPFS) and keeps the record pointers in IndexedDB. A finished run opens read-only in every tab; a running run belongs to the tab that started it through a lock of the Web Locks API, and rename, pin, and delete take the lock and tell the other tabs through a `BroadcastChannel`
  - (+) Run links, the run list, pins, and variants keep their meaning after a reload and in a new tab, with one run semantics in web2 and desktop
  - (-) It adds a mirror with committed file generations and recovery, a "stored, not loaded" state for the run list, quota errors, a storage meter with the `navigator.storage.persist()` state, and cleanup. Recovery must change too: `RunManager::open()` marks every stored running run as interrupted, which is wrong while another tab owns a live run
  - (-) Browsers may clear the store: storage that is not persistent is "best effort", and Safari deletes the script-writable storage of a site after seven days of use without interaction with the site, including the File System API. The run package stays the only durable copy in both models
  - (=) WasmFS with its OPFS back end would also take the files out of tab memory, but it needs Asyncify or JSPI in the single-threaded module and is untested with the wasm-bindgen output
- **Runs that live in the tab** (chosen)
  - (+) The smallest change: the platform fixes that the persistent model needs too (delete route, crash banner, leave warning, missing-run state, run package), and no second storage layer
  - (-) A reload, a closed tab, a crash, or a mobile browser that ends the tab loses every run; `beforeunload` does not fire reliably on mobile platforms and never on a crash
  - (-) The shared interface has two run semantics: a durable history on desktop and the list of one tab on web2

## Reason

The run package is the durable copy in both models, because browsers may evict any storage that the page writes. Persistent browser storage would add a storage layer, cross-tab ownership, and eviction handling to protect a copy that the browser may still delete, while the run package already recreates a run in full on every app. The evidence and the sources are in [kb/reports/2026-11-12-app-ui-redesign-and-web2-audit.md](../reports/2026-11-12-app-ui-redesign-and-web2-audit.md), section "Run model of web2".

A store that outlives the page stays possible later: the record and the package import already recreate a run from files, which a mirror in the origin private file system would reuse.

## Implementation

- `packages/app-ui/src/preferences/storage.ts`: the browser storage keys
- `packages/app-backend/src/runs/package.rs`, `package_import.rs`: export and import of the run package; routes `GET /api/runs/{id}/package` and `POST /api/runs/import` in `packages/app-api/src/routes.rs`
- `packages/app-backend/src/runs/store.rs`, `input_store.rs`: the run folder layout and the inputs stored by content; `DELETE /api/runs/{id}` frees the memory of a run
- `packages/app-ui/src/start/useImportRun.ts`: the import from the start screen and the missing-run state

## Issues

- [kb/issues/M-web2-wasm-panic-loses-the-back-end-state.md](../issues/M-web2-wasm-panic-loses-the-back-end-state.md)
- [kb/issues/M-app-runs-folder-grows-without-retention-or-quota.md](../issues/M-app-runs-folder-grows-without-retention-or-quota.md)
