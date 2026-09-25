# The desktop dev server points Electron at a directory that does not exist

`packages/app-desktop/vite.config.ts` sets `ELECTRON_OVERRIDE_DIST_PATH` to `<checkout>/node_modules/electron/dist` unless the variable is set. With it set, the `electron` package returns `<override>/electron` as the executable path and skips its own download.

## Evidence

- Bun installs the workspace with the isolated linker: the `electron` package lives in `node_modules/.bun/electron@<version>/node_modules/electron` and is linked into `packages/app-desktop/node_modules/electron`; `<checkout>/node_modules/electron` does not exist
- Bun does not run the `postinstall` script of `electron`, because the root `package.json` lists no `trustedDependencies`, so `dist/` and `path.txt` are missing after `bun install`
- `packages/app-desktop/node_modules/electron/index.js` joins the override path with the executable name without checking that it exists

## Impact

- `just desktop` in a fresh checkout starts no Electron binary unless the developer sets `ELECTRON_OVERRIDE_DIST_PATH` to an installed Electron or runs `node packages/app-desktop/node_modules/electron/install.js` first

## Required change

Install the Electron binary as part of the documented setup (for example by trusting the `electron` postinstall script, whose download is checked against the `checksums.json` of the pinned package) and point the override at the linked package, or drop the override so the `electron` package resolves its own `dist/`.

## Validation

- In a fresh worktree, `bun install` followed by `just desktop` (with a display) starts the app without extra environment variables
