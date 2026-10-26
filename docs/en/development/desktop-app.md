---
title: Desktop app
description: Run and package the Electron desktop app, its folders, settings, and example data.
---

`just desktop` starts the Electron app in development mode, with the Vite server of its renderer on the desktop port of the checkout (14007 in the main checkout, `TREETIME_DESKTOP_PORT` overrides it). `just desktop-prod` builds the desktop app for production (the Node addon in the `dist` profile and the bundled renderer) and starts that build from the checkout. In the container both need the host display: `TREETIME_DOCKER_X11=1 ./dev/docker/run just desktop`.

The Node addon (`packages/app-napi`, a Rust `cdylib`) is built by cargo, like the CLI, and copied to `packages/app-napi/app-napi.node`, the `main` file of the package: `just build-napi` builds the `dev` profile for `just desktop`, and `just build-napi dist` the `dist` profile for `just build-desktop` and `just desktop-prod`. `just gen napi-types` generates its TypeScript types (`index.d.ts`).

## Folders and settings

The desktop app keeps all its data in one app folder: `~/.config/treetime` on Linux (`$XDG_CONFIG_HOME/treetime` when set), `~/Library/Application Support/treetime` on macOS, and `%LOCALAPPDATA%\treetime` on Windows. `TREETIME_APP_DIR` names another app folder. The Rust core resolves the folders (`AppPaths` in `packages/app-backend/src/app_paths.rs`), so a later local web mode can use the same ones.

| Content                                              | Default                      | Environment variable      | `settings.yaml` key |
| ---------------------------------------------------- | ---------------------------- | ------------------------- | ------------------- |
| Settings                                             | `<app folder>/settings.yaml` | moves with the app folder | none                |
| Browser profile: cache, local storage, crash reports | `<app folder>/profile`       | `TREETIME_PROFILE_DIR`    | `paths.profile`     |
| Runs                                                 | `<app folder>/runs`          | `TREETIME_RUNS_DIR`       | `paths.runs`        |
| Logs and crash diagnostics                           | `<app folder>/logs`          | `TREETIME_LOGS_DIR`       | `paths.logs`        |
| Example datasets and configurations                  | `<app folder>/examples`      | `TREETIME_EXAMPLES_DIR`   | `paths.examples`    |

The `--app-dir=<dir>` switch of the app names the app folder for one start, relative to the folder the app starts from: `treetime-desktop --app-dir=/tmp/treetime-check`. It takes precedence over `TREETIME_APP_DIR`, which takes precedence over the development folders below and the platform default. Use it, not Chromium's `--user-data-dir`, which moves only the browser profile, so two processes would still share `settings.yaml` and `runs/`. One instance of the app runs per browser profile folder: a second start with the same profile folder focuses the window of the first one. `TREETIME_PROFILE_DIR` or `paths.profile` can make two app folders share one profile folder, and then only one of them runs at a time.

An environment variable takes precedence over the settings key, and the key over the default. A relative path in `settings.yaml` is relative to the app folder; a relative path in a variable is relative to the working directory. `just desktop` and `just desktop-prod` use `tmp/app/treetime-dev` and `tmp/app/treetime-prod` of the checkout as app folder, so development never touches an installed app, and `data/` of the checkout as examples folder; `TREETIME_DESKTOP_DEV_DIR`, `TREETIME_DESKTOP_PROD_DIR`, and the variables above in `.env` choose other folders (see `.env.example`).

## Example data

Every app takes its example data from one index, `data/index.json` (`DataIndex` in `packages/app-backend/src/datasets.rs`): every file that git tracks in `data/` except hidden files and folders, the datasets with their input files, the example configurations with the files they read, and the pipelines. The examples panel, the command palette, and the "Example files" action of each file input read the index.

- **Web apps**: the web app and web2 serve `data/index.json` and the files it lists from their own origin, through the Vite plugin `exampleData()` (`packages/app-ui/build/example-data-vite.ts`), which serves them in the dev server and copies them into `dist/data/` in the build. A chosen example file is uploaded into the run like a user file, and an example configuration is loaded with the paths of its uploaded files
- **Desktop app**: an installed app starts with an empty examples folder. The example dialogs and the command palette then offer "Download examples": the back end downloads `examples.zip` of the app's release from `neherlab/treetime-nightly` and unpacks it into the examples folder in the background (`POST /api/examples/download`, with the progress on the app event stream), so a reload of the window does not stop it. The archive holds exactly the files of `data/index.json` and the index itself, and `GET /api/examples` returns the examples folder with its index. A nightly build downloads the archive of its own release; a development or release build downloads the archive of the latest release. The app reads the example files where they are, without an upload

No CLI command downloads examples; CLI users download `examples.zip` from the release.

`settings.yaml` holds the folders (`paths`) and the preferences of the user interface (`ui`: theme, sidebar width, and the unfinished form):

```yaml
paths:
  runs: "/data/treetime-runs"
ui:
  theme: "dark"
```

Every key is optional, and an unknown key is an error that names the file. When `settings.json` exists instead, the app reads and writes JSON; both files at once is an error. The app reads the folders when it starts. The user changes the runs folder in the app, with the folder button in the header or "Change runs folder" in the command palette; the back end then restarts to open it, which interrupts runs that are computing. A runs folder set by `TREETIME_RUNS_DIR` cannot be changed in the app. Only the desktop app serves the settings routes (`/api/app-settings`, `/api/workspace`): `treetime-server` leaves them out, so a shared server never mixes the preferences of its visitors.

## Packages

`just package-desktop [target]` packages the desktop app with [electron-builder](https://www.electron.build) into `.out/treetime-desktop-<target>.<ext>`:

- Linux (`x86_64-unknown-linux-gnu`, `aarch64-unknown-linux-gnu`): an AppImage, which runs without installation
- Windows (`x86_64-pc-windows-gnullvm`, `aarch64-pc-windows-gnullvm`): a zip archive; users extract it and start `TreeTime.exe`
- macOS (`x86_64-apple-darwin`, `aarch64-apple-darwin`): a dmg

The recipe builds the addon for the native target `x86_64-unknown-linux-gnu` itself; the other targets take the addon from `just cross-napi` (see [Cross-compilation and releases](/development/releases)). Linux and Windows packages build in the container; electron-builder downloads Electron, the AppImage tools, and 7-Zip into `.build/desktop/cache/`. macOS packages build only on macOS, because the dmg tools and `codesign` exist only there. `dev/desktop/package --help` describes the steps: `dev/desktop/stage` puts the bundles, a minimal `package.json`, and the addon into `.build/desktop/<target>/app`, electron-builder packages that directory with `packages/app-desktop/electron-builder.config.ts`, and `dev/desktop/check` checks the Electron fuses and that the addon sits outside `app.asar`. `dev/desktop/start-test <target>` starts a package and fails when it exits early or logs an error.

The packages set the Electron [fuses](https://www.electronjs.org/docs/latest/tutorial/fuses) that turn off running the app as plain Node.js (`RunAsNode`, `NODE_OPTIONS`, `--inspect`) and load the app only from its integrity-checked `app.asar`. They are not signed with a developer certificate: on macOS they carry an ad-hoc signature, and macOS (Gatekeeper) and Windows (SmartScreen) ask the user to confirm the first start. The app has the default Electron icon.

On Ubuntu 24.04 and newer, AppArmor blocks the unprivileged user namespaces that the Chromium sandbox needs. The AppImage cannot install an AppArmor profile, so its start script detects the block and starts the app with `--no-sandbox`: on those systems the app runs without the Chromium sandbox. The AppImage uses the static AppImage runtime, which needs only `fusermount3` (installed by default on Ubuntu) and not `libfuse2`.
