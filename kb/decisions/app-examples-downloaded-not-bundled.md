# Example datasets are downloaded, not bundled

The example datasets in `data/` ship as one `examples.zip` asset of each GitHub release in `neherlab/treetime-nightly`. Installed desktop apps and CLI users download that archive when they ask for it. The desktop packages and the CLI binaries contain no datasets.

## Behavior

- **Release asset**: the nightly publish job builds `examples.zip` from `data/`. Archive entries start at the dataset folders (`zika/20/tree.nwk`), without the smoke fixtures in `data/smoke/` and without the development scripts in `data/`
- **Release choice**: a nightly build downloads the archive of its own release tag. Development and release builds download the archive of the latest release
- **CLI**: `treetime examples get --output-dir <dir>` downloads and unpacks the archive into an empty or missing folder. `--url` downloads another archive
- **Desktop**: when the examples folder holds no datasets and no example configs, the example dialogs and the command palette offer a "Download examples" button. The back end downloads and unpacks the archive into the examples folder and reports progress on the app event stream. The download runs in the back end, so it continues when the window reloads
- **Web app**: the hosted web app copies `data/` into its image and does not offer the download
- **Offline start**: an installed app starts without the datasets. The rule to bundle every runtime asset (fonts, scripts, styles) does not cover example datasets, because they are optional data that no feature needs to work

## Reason

Bundling would put a catalog of about 22 MB into every desktop package and every CLI binary, while most users run their own data.

## Implementation

- `.github/workflows/nightly.yml`, `dev/publish-nightly`: the release asset
- `packages/app-commands/src/examples_download.rs`: the downloader that the CLI and the desktop back end share
- `packages/app-cli/src/cli/examples.rs`: `treetime examples get`
- `packages/app-server/src/app_settings_routes.rs`: `POST /api/examples/download` and `GET /api/examples/download`, local apps only
