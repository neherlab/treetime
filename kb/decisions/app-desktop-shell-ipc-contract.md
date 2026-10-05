# Desktop shell messages and the save trust boundary

The Electron renderer talks to two processes. It reaches Rust over HTTP through a message port to the utility process, with the generated client. It reaches the Electron main process over IPC for native functions: file and folder dialogs, the native theme, restarting the back end, and saving run files.

## Message types

- **Values that reach Rust or carry back-end data** are Rust types that derive their JSON schema. TypeScript gets their types and zod validators only from `just gen`, and every transport validates incoming data with the generated zod
- **Messages between the renderer and the main process that never reach Rust** use hand-written zod schemas. They are declared once, in the channel table of `packages/app-ui/src/host.ts`, and the main process checks every request against its schema before the handler runs. Any part of such a message that reaches Rust or carries back-end data is built from the generated schemas: the `save-run` request takes its run id and file path from the generated `runsSave` schemas, the theme message is the generated `UiTheme`, and the `backend-stopped` event carries the generated `ErrorResponse`
- **One interface**: `interface Host` derives its method types from the channel table. The desktop preload script exposes one host object; the web app has no host. `packages/app-ui` receives the host through one React context

## Save trust boundary

Saving a run file writes to a path that the user picks, so the renderer never names a file-system destination:

- The renderer sends `save-run` with the run id, the optional output file path inside the run, and a suggested file name
- The main process shows the save dialog, then sends `POST /api/runs/{id}/save` with the destination to the back end itself, over a message port that only the main process holds
- The back end builds two routers from one API: the host router, with the save route, and the renderer router, without it. The utility process learns the scope of each port from the control message of the main process, which only the main process can post. A renderer request to the save route answers 404
- In the web app, the same buttons start a browser download of the file or the archive

## Reason

One declaration per message removes the copies of the dialog and save types that lived in several packages, and the save boundary keeps a compromised renderer from writing to arbitrary files through the back end.

## Implementation

- `packages/app-ui/src/host.ts`: channel table and `interface Host`
- `packages/app-desktop/src/ipc-main.ts`, `packages/app-desktop/src/ipc-renderer.ts`: the only users of `ipcMain` and `ipcRenderer`
- `packages/app-napi/src/backend.rs`: host and renderer routers
- `packages/app-server/src/app_settings_routes.rs`: `runsSave`
