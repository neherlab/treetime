# The desktop dataset catalog depends on the working directory

The desktop back end reads the example datasets and example configurations from `DATA_DIR`, or from `data` relative to the working directory of the process when `DATA_DIR` is unset (`DesktopService::open` in `packages/app-napi/src/backend.rs`). `just desktop` and `just desktop-prod` change into the checkout first (`TREETIME_PROJECT_ROOT` in `packages/app-desktop/src/main.ts`), so they find `data/` of the checkout. An installed app starts in whatever directory the system or the user launched it from, and the packages do not bundle `data/`.

## Evidence

- `DATA_DIR_ENV` and `DEFAULT_DATA_DIR = "data"` in `packages/app-napi/src/backend.rs`; the path is used as given, without resolving it against the app or the user folders
- `discover_datasets` in `packages/app-datasets/src/lib.rs` treats a missing directory as an empty catalog, so the app shows no examples instead of an error
- `dev/desktop/stage` copies only the Vite bundles and the Node addon into the package; `packages/app-desktop/electron-builder.config.ts` adds no extra resources

## Impact

- An installed app usually shows an empty example catalog
- An app started from a directory that happens to contain a `data` folder lists that folder's files as examples

## Options

- Bundle a curated set of example datasets as extra resources of the package and resolve `DATA_DIR` against the resources folder of the app
- Resolve the catalog folder in `AppPaths` (`packages/app-commands/src/app_paths.rs`), next to the other folders of the app, and let the user add datasets there
- Download example datasets on demand into the data folder of the app

## Locations

- `packages/app-napi/src/backend.rs`: `DesktopService::open`
- `packages/app-desktop/src/main.ts`: working directory of the main and back-end processes
- `dev/desktop/stage`, `packages/app-desktop/electron-builder.config.ts`: package contents
