# The desktop app does not set Electron fuses

Electron fuses are build-time switches in the Electron binary that turn off features an attacker can use to run code with the application's privileges, for example `RunAsNode` (the `ELECTRON_RUN_AS_NODE` environment variable turns the app into a plain Node.js interpreter) and `EnableEmbeddedAsarIntegrityValidation`. The Electron security checklist recommends flipping them for every shipped app ([security tutorial, item 19](https://www.electronjs.org/docs/latest/tutorial/security), [fuses](https://www.electronjs.org/docs/latest/tutorial/fuses)).

## Evidence

- The repository has no packaging step for the desktop app: no electron-builder, Electron Forge, or `@electron/fuses` configuration, and `just build-desktop` produces only the Vite bundles and the addon
- Fuses are written into the packaged Electron binary, so they can only be set by a packaging step; the development binary under `node_modules` must keep its defaults

## Impact

- A packaged build made by hand keeps `RunAsNode` and the other fuses at their permissive defaults

## Required change

When a packaging step is added, flip the fuses in it: turn off `RunAsNode`, `EnableNodeOptionsEnvironmentVariable` and `EnableNodeCliInspectArguments`, and turn on `EnableCookieEncryption`, `EnableEmbeddedAsarIntegrityValidation` and `OnlyLoadAppFromAsar`. Verify that the back end still starts in its utility process with `RunAsNode` off, because the Electron documentation does not state whether `utilityProcess` depends on it.

## Validation

- `npx @electron/fuses read --app <packaged binary>` lists the fuses as set
- The packaged app starts, and the back end answers in its utility process
