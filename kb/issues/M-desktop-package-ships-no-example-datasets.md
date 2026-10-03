# The desktop package ships no example datasets

The desktop app lists the example datasets and configurations of its examples folder: `<app folder>/examples`, `paths.examples` in `settings.yaml`, or `TREETIME_EXAMPLES_DIR` (`AppPaths` in `packages/app-commands/src/app_paths.rs`). `just desktop` and `just desktop-prod` point it at `data/` of the checkout. An installed app starts with an empty examples folder, because the packages do not bundle any datasets, so its example list is empty until the user copies datasets there.

## Evidence

- `DesktopService::open` in `packages/app-napi/src/backend.rs` creates the examples folder and passes it to the dataset catalog
- `discover_datasets` in `packages/app-datasets/src/lib.rs` lists the datasets below that folder
- `dev/desktop/stage` copies only the Vite bundles and the Node addon into the package; `packages/app-desktop/electron-builder.config.ts` adds no extra resources

## Impact

- New users of the desktop app see no examples to start from; the web app shows the examples of the server's `--data-dir`

## Options

- Bundle a curated set of small datasets as extra resources of the package, and copy them into the examples folder at the first start
- Download example datasets on demand into the examples folder
- List bundled read-only examples next to the examples folder of the user

## Locations

- `packages/app-napi/src/backend.rs`: `DesktopService::open`
- `dev/desktop/stage`, `packages/app-desktop/electron-builder.config.ts`: package contents
