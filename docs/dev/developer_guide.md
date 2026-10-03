# Developer guide

This guide describes how to set up a development environment, build and test TreeTime, and maintain the project. It assumes basic familiarity with TreeTime from a user perspective.

The `justfile` is the entry point for every routine task. `just` lists the recipes by group, and each recipe has a one-line description.

## Setup

TreeTime builds in two ways: in the build container (recommended), or directly on the host. Both run the same `just` recipes with the same tool versions, which `.config/mise.toml` pins and `.config/mise.lock` locks by URL and checksum.

### Build container (recommended)

Requirements: Docker with buildx, git, and bash (the stock bash 3.2 of macOS works).

```bash
git clone https://github.com/neherlab/treetime
cd treetime
./dev/docker/run just setup
./dev/docker/run just check
```

`./dev/docker/run <command>` runs a command in the container, and without a command it opens a shell. The first run builds the image, which takes a while; later runs reuse it until one of its build inputs changes. Notes:

- The checkout is mounted at its host path, and the git metadata is mounted read-only. Commit on the host
- Build output goes to `.build/container/`, the cargo home and caches to `.cache/`
- Commands run as your user, with all Linux capabilities dropped
- On Apple Silicon and other arm64 hosts the image runs under x86_64 emulation, which is slower
- The container uses host networking, so the dev servers are reachable from the host browser. On Docker Desktop, turn on host networking in the settings

### Host

Requirements: [mise](https://mise.jdx.dev), [rustup](https://rustup.rs), bash 4 or later (on macOS: `brew install bash`), a C toolchain with gfortran, OpenBLAS, and libclang. On Debian and Ubuntu:

```bash
sudo apt-get install build-essential gfortran libclang-dev libfontconfig-dev libopenblas-dev pkg-config
mise install
just setup
just check
```

Host builds go to `.build/host/`. Host builds of the Rust code work on Linux x86_64 only; on macOS, build in the container. Every build links OpenBLAS statically: `.cargo/config.toml` makes the linker take `libopenblas.a` even where `libopenblas-dev` also installs the shared library, so no binary and no Node addon depends on a system OpenBLAS. On macOS, the pinned `openblas-src` links Homebrew's OpenBLAS, which also needs the gfortran and OpenMP runtimes, and the configuration passes no link flags for them. The dylint and IQ-TREE tools are available on Linux only; run the recipes that need them in the container on other hosts.

### Machine settings

Optional settings go into the gitignored `.env` in the checkout; `.env.example` lists them:

- `KACHE_STORE`: directory of the [kache](https://github.com/kunobi-ninja/kache) compiler cache stores. Every build and clippy pass then compiles through kache, which shares compiled crates across the worktrees of this project. Workspace crates of incremental builds (`dev`, tests, `release`, clippy) keep their incremental state, so an edit rebuilds as fast as without kache; kache serves their dependencies and every build without incremental state (CI, `dist`, `profiling`, `bench`, cross builds). Dylint, hawk, and coverage compile without kache. The directory holds one store per kache version, environment (`host` or `docker`), and pass (`build`, `clippy`, or `cross-<target>`)
- `KACHE_MAX_SIZE`: size limit of each kache store, 100 GiB when unset
- `SMOKE_STORE`: absolute directory of a smoke snapshot store shared by all checkouts of this project on the machine, so each baseline is built once (see [Snapshot store](#snapshot-store)). Unset: each checkout keeps its own `snapshots/`
- `TREETIME_PORTLESS`: use of the [portless](https://github.com/vercel-labs/portless) proxy by the web dev server (see [Web app](#web-app))

## Everyday commands

| Task                                             | Command                              |
| ------------------------------------------------ | ------------------------------------ |
| List the recipes                                 | `just`                               |
| Fast checks (format, clippy, TypeScript)         | `just check`                         |
| Every check; must pass before merging            | `just check-all`                     |
| Fast lint fixes and format (stage changes first) | `just fix`                           |
| Every lint fix, dylint included, and format      | `just fix-all`                       |
| Build                                            | `just build` (`just b`)              |
| Run the CLI                                      | `just run treetime ancestral --help` |
| Rust tests, optionally filtered                  | `just test-rs [filter]` (`just t`)   |
| TypeScript tests                                 | `just test-ts`                       |
| Clippy                                           | `just lint-rs` (`just l`)            |
| TypeScript lints                                 | `just lint-ts`                       |
| Format one toolchain                             | `just fmt-rs`, `fmt-ts`, `fmt-other` |

In the container, prefix each command with `./dev/docker/run`.

Recipe names follow one scheme:

- **Leaf recipes** run one tool in one mode. The suffix `-rs` or `-ts` names the toolchain, Rust or TypeScript; tools that serve one toolchain only keep their own name, such as `dylint`, `hawk`, or `knip`
- **Combined recipes** without a language suffix, such as `lint`, `test`, `fmt`, or `fix`, run the leaves of both toolchains and call no tool themselves
- **The suffix `-all`** adds the slow tools to the fast set of the same name: `fix-all` adds the dylint fixes, `lint-all` the custom lint libraries, hawk, and the dependency and config lints, `test-all` the tests of the lint libraries, and `check-all` the full gate

While working on one toolchain, run its leaves (`just fmt-rs`, `just l`, `just t`) and leave the slow tools to `check-all`.

`check` and `check-all` run their checks in parallel through `dev/run-checks`, keep going past failures, and list the failed checks at the end. Each check writes its output to `tmp/checks/<check>.log`. Warnings fail the checks, and a missing tool is a failure, not a skipped check. The full gate consists of groups, one CI job each; `just check-group <group>` runs one group serially, as its CI job does.

### Examples

Replace `$v` with a dataset path under `data/`, for example `flu/h3n2/20`:

```bash
just run treetime ancestral --method-anc=marginal --tree=data/$v/tree.nwk --alignment=data/$v/aln.fasta.xz --output-all=tmp/ancestral/$v
just run treetime clock --tree=data/$v/tree.nwk --metadata=data/$v/metadata.tsv --output-all=tmp/clock/$v
just run treetime timetree --tree=data/$v/tree.nwk --metadata=data/$v/metadata.tsv --alignment=data/$v/aln.fasta.xz --output-all=tmp/timetree/$v
```

## Project structure

```
packages/
  app-cli/           CLI binary (treetime)
  app-commands/      Command configs and runners shared by the CLI, the server, and the addon
  app-contracts/     TypeScript types and schemas of the bridge between the apps and the Rust code
  app-datasets/      Bundled example datasets
  app-desktop/       Electron shell of the desktop app
  app-napi/          Node addon of the desktop app (Rust, napi-rs)
  app-output/        Output writers shared by the CLI, the server, and the addon
  app-server/        HTTP API server of the web app (treetime-server)
  app-ui/            React components shared by the desktop and web apps
  app-web/           Web app (Vite)
  legacy/            TreeTime v0, the Python reference implementation
  schemas/           Generated JSON schemas of the CLI input
  treetime/          Core library
  treetime-*/        Supporting crates of the core library
  util-*/            File format libraries
dev/                 Development scripts and the container setup
test_scripts/        Python research scripts and notebooks
kb/                  Knowledge base
data/                Example datasets
```

## Apps

The TypeScript apps (`app-ui`, `app-web`, `app-desktop`) are a user interface around the Rust core. Algorithms, file parsers, readers and writers, configuration handling, and domain rules live in Rust and are tested there. The apps render what the Rust operations return and keep only presentation logic: layout, interaction state, and formatting for display. One implementation of each rule keeps the CLI, the web app, and the desktop app from diverging, and one test suite covers it.

- **CLI knowledge**: flags, value syntax, setting groups, defaults, the equivalent command line, and the YAML config come from clap and the config types. The setting catalog in the OpenAPI document (`x-setting-catalog`) describes every setting, and `run-config` and `check-config` render the command line and YAML
- **Checks**: `check-config` classifies configuration problems and input facts into checks that block the run, warn, or advise, with the settings that fix them
- **Types**: every value that crosses from Rust to TypeScript is a Rust type with a derived schema in the OpenAPI document; the apps read it through the generated types and zod schemas in `app-contracts`, never through hand-written ones

### Interface components

The interface is built from [shadcn/ui](https://ui.shadcn.com) components in the `base-vega` style, which run on Base UI. The components are vendored into `packages/app-ui/src/ui/`, one file per component, and the app composes them; it does not style raw elements. TreeTime owns the vendored files: they pass the same lint rules as the rest of the source and keep only the parts the app uses. `packages/app-ui/src/ui/UPSTREAM.md` records their source, the rewrites that turn a registry file into a component, and the local changes to re-apply when a component is refreshed.

The theme tokens (colors, radius, fonts) live in `packages/app-ui/src/theme.css`: green-grey surfaces with a deep teal accent, Lato for text and IBM Plex Mono for code. `packages/app-ui/src/ui/shadcn.css` is the stylesheet of the `shadcn` npm package that defines the variants the components use.

Icons come from [Iconify](https://icon-sets.iconify.design) sets. `unplugin-icons` compiles each imported icon into an inline SVG React component at build time, so the bundle holds only the icons in use and the apps load nothing from an icon host at runtime. An icon is a default import from the virtual `~icons/<set>/<name>` module:

```tsx
import XIcon from "~icons/lucide/x";

<XIcon aria-hidden className="size-4" />;
```

Lucide is the base set and matches the stroke style of the shadcn components. When Lucide has no fitting icon, add the Iconify set that has one: install its `@iconify-json/<set>` package as a pinned devDependency of `packages/app-ui` and import from `~icons/<set>/<name>`; the plugin configuration in `packages/app-ui/build/icons-vite.ts` needs no change. Icons render at 24 by 24 pixels unless a class or a parent sets their size. A decorative icon takes `aria-hidden`, and an icon that carries meaning on its own takes an `aria-label`. `lucide-react` is banned by the lint rules because its components duplicate the compiled icons.

Library hooks cover the behavior around the components: `@tanstack/react-hotkeys` for keyboard shortcuts and their platform labels, `react-dropzone` for file drops and the file dialog, `use-stick-to-bottom` for the following log, `@tanstack/react-table` for sortable tables, `@tanstack/react-pacer` for debouncing, Base UI Autocomplete for the command palette, and `@mantine/hooks` for the clipboard, the size of the Auspice panel, and the mobile media query.

### Web app

`just up` starts the API server and the Vite dev server in the foreground until Ctrl-C. `just health` and `just status` report whether they run and which commit they were started from.

`just serve` builds the web app and the API server for production and runs the server in the foreground until Ctrl-C. The server serves the built web app (`packages/app-web/dist`, through the `STATIC_DIR` environment variable) and the API from one port.

Each checkout resolves its own ports, so worktrees run side by side; `TREETIME_API_PORT`, `TREETIME_WEB_PORT`, and `TREETIME_SERVE_PORT` override them.

#### Storage

The web app has no database. The server keeps each run in its own folder `<runs-dir>/<run-id>/`: the run record `run.json`, the event log `events.jsonl`, the uploaded files in `inputs/`, and the results in `out/`. The browser keeps the preferences of the user interface (theme, sidebar width, and the unfinished form) in `localStorage`, under the key `treetime-preferences`. The server reads example datasets from `--data-dir`.

Each app mode keeps its runs in its own directory:

- `tmp/app/web-dev/runs`: development servers (`just up`)
- `tmp/app/web-prod/runs`: production server (`just serve`)
- `tmp/app/treetime-dev/runs`: desktop app in development mode (`just desktop`); see [Desktop app](#desktop-app) for its other folders
- `tmp/app/treetime-prod/runs`: production build of the desktop app started from the checkout (`just desktop-prod`)
- `~/.local/share/treetime/runs` on Linux, `~/Library/Application Support/treetime/runs` on macOS, `%LOCALAPPDATA%\treetime\runs` on Windows: installed desktop app, unless the user chose another runs folder

A deployment passes `--data-dir` and `--runs-dir` to `treetime-server` itself.

When the portless proxy runs on the machine, the web server also registers `https://treetime.localhost` in the main checkout and `https://<branch>.treetime.localhost` in a linked worktree.

#### Event streams and HTTP/2

The web app follows changes through server-sent event streams: `GET /api/events` for changes to any run, and `GET /api/runs/{id}/events` for each run it follows. Every open stream holds one connection. Over HTTP/1.1, browsers allow at most 6 connections per host, shared by all tabs, so more open streams stall every further request to the server. Over HTTP/2, all requests share one connection.

The API server speaks HTTP/1.1 without TLS, and browsers use HTTP/2 only over TLS. A deployment therefore puts a reverse proxy in front of the server that terminates TLS and serves HTTP/2 to the browser, for example nginx with `listen 443 ssl; http2 on;`. In development, the portless proxy (`https://treetime.localhost`) does this by default. Without such a proxy, at most 6 streams stay open per host across all tabs.

#### Caching, compression, and allowed hosts

`treetime-server` serves the UI build (`STATIC_DIR`) and the API with these rules, so a proxy in front of it needs no caching or compression settings of its own:

- **UI files**: the file names in `assets/` contain a hash of their content, so they are cached for a year (`Cache-Control: public, max-age=31536000, immutable`). Every other file, `index.html` included, is revalidated on each use (`no-cache`), so a new deploy reaches the browser at the next page load. Only page navigations (`Accept: text/html`) fall back to `index.html`; a missing asset or API path answers 404
- **Compression of UI files**: `bun run build:web` writes a gzip (level 9) and a brotli (quality 9) copy next to each text file of the build, and the server sends the copy the browser accepts, with `Vary: Accept-Encoding`. Brotli quality 11 would make the files about 8 % smaller but adds seconds to every build
- **Compression of API responses**: zstd, or gzip for browsers without zstd, at the fastest level. Run outputs of large alignments shrink far more with zstd than with gzip, because gzip only finds repeats within 32 KB, less than one genome. Event streams and already compressed files (`.xz`, `.gz`, `.zst`, fonts) stay uncompressed
- **API caching**: API responses are not stored (`no-store`). The results, the Auspice document, and the comparison of finished runs carry an `ETag`, so a browser that asks again gets `304 Not Modified` without the server rebuilding the answer. The tag changes when the server restarts
- **Downloads**: output files and run archives stream to the browser, so the server holds only a small buffer per download. Archives store their files uncompressed and rely on the HTTP compression above
- **Allowed hosts**: the server answers only requests for `localhost`, `*.localhost`, `127.0.0.1`, `[::1]`, and the hosts given with `--allowed-host` (repeatable) or `ALLOWED_HOSTS` (comma-separated); other hosts get HTTP 403. This stops other websites from reaching a local server through DNS rebinding. The server sends no CORS headers, because the UI calls it from its own origin

### Desktop app

`just desktop` starts the Electron app in development mode. `just desktop-prod` builds the desktop app for production (the Node addon in the `dist` profile and the bundled renderer) and starts that build from the checkout. In the container both need the host display: `TREETIME_DOCKER_X11=1 ./dev/docker/run just desktop`.

The Node addon (`packages/app-napi`, a Rust `cdylib`) is built by cargo, like the CLI, and copied to `packages/app-napi/app-napi.node`, the `main` file of the package: `just build-napi` builds the `dev` profile for `just desktop`, and `just build-napi dist` the `dist` profile for `just build-desktop` and `just desktop-prod`. `just gen napi-types` generates its TypeScript types (`index.d.ts`).

#### Folders and settings

The Rust core resolves the folders of the desktop app (`AppPaths` in `packages/app-commands/src/app_paths.rs`), so a later local web mode can use the same ones. Each folder is named `treetime`:

| Content                                         | Linux                          | macOS                                         | Windows                         |
| ----------------------------------------------- | ------------------------------ | --------------------------------------------- | ------------------------------- |
| Chromium profile and `settings.yaml`            | `~/.config/treetime`           | `~/Library/Application Support/treetime`      | `%LOCALAPPDATA%\treetime`       |
| Runs, unless the settings name another folder   | `~/.local/share/treetime/runs` | `~/Library/Application Support/treetime/runs` | `%LOCALAPPDATA%\treetime\runs`  |
| Logs and crash diagnostics                      | `~/.local/state/treetime/logs` | `~/Library/Logs/treetime`                     | `%LOCALAPPDATA%\treetime\logs`  |

On Linux, `XDG_CONFIG_HOME`, `XDG_DATA_HOME`, and `XDG_STATE_HOME` move these folders. `TREETIME_APP_DIR` replaces all of them with one folder that holds `profile/`, `settings.yaml`, `runs/`, and `logs/`. `just desktop` and `just desktop-prod` set it to `tmp/app/treetime-dev` and `tmp/app/treetime-prod` of the checkout, so development never touches the profile of an installed app; `TREETIME_DESKTOP_DEV_DIR` and `TREETIME_DESKTOP_PROD_DIR` in `.env` choose other folders (see `.env.example`).

`settings.yaml` holds the runs folder (`workspace`) and the preferences of the user interface (`ui`: theme, sidebar width, and the unfinished form). Every key is optional, and an unknown key is an error that names the file. When `settings.json` exists instead, the app reads and writes JSON; both files at once is an error. The user changes the runs folder in the app, with the folder button in the header or "Change runs folder" in the command palette; the back end then restarts to open it, which interrupts runs that are computing. Only the desktop app serves the settings routes (`/api/app-settings`, `/api/workspace`): `treetime-server` leaves them out, so a shared server never mixes the preferences of its visitors.

#### Packages

`just package-desktop [target]` packages the desktop app with [electron-builder](https://www.electron.build) into `.out/treetime-desktop-<target>.<ext>`:

- Linux (`x86_64-unknown-linux-gnu`, `aarch64-unknown-linux-gnu`): an AppImage, which runs without installation
- Windows (`x86_64-pc-windows-gnu`): a zip archive; users extract it and start `TreeTime.exe`
- macOS (`x86_64-apple-darwin`, `aarch64-apple-darwin`): a dmg

The recipe builds the addon for the native target `x86_64-unknown-linux-gnu` itself; the other targets take the addon from `just cross-napi` (see [Cross-compilation and releases](#cross-compilation-and-releases)). Linux and Windows packages build in the container; electron-builder downloads Electron, the AppImage tools, and 7-Zip into `.build/desktop/cache/`. macOS packages build only on macOS, because the dmg tools and `codesign` exist only there. `dev/desktop/package --help` describes the steps: `dev/desktop/stage` puts the bundles, a minimal `package.json`, and the addon into `.build/desktop/<target>/app`, electron-builder packages that directory with `packages/app-desktop/electron-builder.config.ts`, and `dev/desktop/check` checks the Electron fuses and that the addon sits outside `app.asar`. `dev/desktop/start-test <target>` starts a package and fails when it exits early or logs an error.

The packages set the Electron [fuses](https://www.electronjs.org/docs/latest/tutorial/fuses) that turn off running the app as plain Node.js (`RunAsNode`, `NODE_OPTIONS`, `--inspect`) and load the app only from its integrity-checked `app.asar`. They are not signed with a developer certificate: on macOS they carry an ad-hoc signature, and macOS (Gatekeeper) and Windows (SmartScreen) ask the user to confirm the first start. The app has the default Electron icon.

On Ubuntu 24.04 and newer, AppArmor blocks the unprivileged user namespaces that the Chromium sandbox needs. The AppImage cannot install an AppArmor profile, so its start script detects the block and starts the app with `--no-sandbox`: on those systems the app runs without the Chromium sandbox. The AppImage uses the static AppImage runtime, which needs only `fusermount3` (installed by default on Ubuntu) and not `libfuse2`.

## Generated files

The JSON schemas, the OpenAPI document, its TypeScript client, and the CLI reference documentation are generated and committed. `just gen` regenerates them, and `just generated-check`, part of `check-all`, fails when a committed copy is stale. Never edit them by hand.

`just fixtures-check` checks the reference and golden-master test fixtures against `dev/registry/reference-files.toml`.

## Testing against the reference

`dev/smoke` (host, needs Docker) runs the CLI over a matrix of commands, datasets and flag variants, stores the outputs in a snapshot named after the source tree that built the binary, and compares them byte for byte with a baseline snapshot, by default the one of the tip of the base branch `rust`. It reports crashes, timeouts, missing outputs and changed outputs, with the command, stderr and the changed values of each case. The cases are defined in `dev/smoke.toml`, one row per command variant with its flags, datasets and expected failures; `./dev/smoke --overview` prints the cases per command and variant. `./dev/smoke --help` describes every option, the snapshot layout and the exit codes.

- `just smoke` compares the quick tier with `rust`
- `just smoke-all` compares every case with `rust`
- `just smoke-run` runs the cases without a baseline: crash, timeout and output checks only
- `just smoke-failed` runs again the cases that did not pass in the last run
- `just smoke-prune` deletes the dirty snapshots of the checkout other than the current one, and old snapshots that no branch needs

```bash
./dev/smoke                    # quick tier against rust
./dev/smoke --tier full        # every case
./dev/smoke --against <ref>    # another commit, branch or snapshot id as the baseline
./dev/smoke --no-compare       # statuses only
./dev/smoke --only <regex>     # cases whose id matches
./dev/smoke --rerun-failed     # cases that did not pass in the last run
./dev/smoke --list             # print the selected case ids
```

Other options:

- `--rerun`: run the selected cases again even when their stored results are reusable
- `--no-build`: use the existing `.out/treetime` for the current snapshot
- `--check-determinism`: run each case a second time with `treetime --jobs=1` and mark the case nondeterministic when the outputs differ
- `--jobs N`, `--timeout-scale X`: parallel cases, and a factor on every case timeout
- `--memory-budget GIB`: memory the running cases may use together, 75% of the available memory by default. Every treetime process is capped at the budget, so a case that needs more fails with an allocation error instead of exhausting the machine

Tiers:

- `quick` (default): datasets of at most 100 sequences, the help texts and a few larger timetree cases, about 6 minutes
- `full`: every case, about 40 minutes

### Snapshot store

Without settings, each checkout keeps its own store: `snapshots/` at the checkout root (git-ignored), with the temporary baseline worktrees and the lock files under `tmp/smoke/`. With `SMOKE_STORE` set to an absolute directory in `.env` (see [Machine settings](#machine-settings)), all checkouts on the machine share one store with `snapshots/`, `worktrees/` and `locks/` in that directory. A baseline is then built once per machine, and every worktree reuses its binary and case results.

- **Snapshot id**: the first 10 hex digits of the git tree id of the commit, or `<tree>+dirty-<hash>` for a working tree with uncommitted or untracked changes. The id depends on content only, so commits with the same tree, such as a branch and the target it was fast-forwarded into, share one snapshot. `--against` resolves a git ref to the tree of its commit, and also accepts an existing snapshot id
- **Binaries**: the current binary is built in the checkout unless its snapshot already has one. A git ref baseline is built in a temporary worktree of the store with a copy of the checkout's `.env`, so it uses the same compiler cache
- **Reuse**: a stored case is reused when its command, expectation, declared outputs, working directory and input files are unchanged. A snapshot whose binary changed loses its stored cases
- **Locks**: a run holds an exclusive lock on a snapshot while it builds and runs it, and a shared lock while it reads the results for its report. A second run that needs the same snapshot logs the pid and checkout of the holder, waits, and then reuses the stored results. Runs on different snapshots never wait for each other. The lock is released when its holder exits, also on a crash, and an interrupted build leaves no binary behind
- **Pruning**: `--prune` deletes the dirty snapshots that the current checkout wrote, other than its current one, and the snapshots of committed trees older than 7 days whose tree is not the tip of any local branch. It skips snapshots that another run holds, so one worktree never deletes the work of another

Layout of one snapshot:

```text
snapshots/<id>/
  snapshot.json             commit, tree, dirty hash, ref, checkout, binary SHA-256, creation time
  bin/treetime              the binary that produced the snapshot
  last-run.json             statuses of the last run and its checkout (for --rerun-failed)
  cases/<case-id>/
    out/                    raw outputs, never modified
    stdout.txt, stderr.txt
    cmd.sh                  the exact command, runnable from the checkout root in the container
    result.json             status, exit code, wall and CPU time, peak memory, case key
  compare/<baseline-id>/    report.md and report.tsv of a comparison
  report.md, report.tsv     report of a --no-compare run
```

`report.md` groups failing cases by likely cause and lists value-level differences. Before the comparison, only the volatile fields `meta.updated` of `*.auspice.json` and `generated_by.version` of `*.augur-node-data.json` are removed; numeric differences are reported, never tolerated. A `changed` case, a `regressed` case, a `still-failing` case without a declared expectation, and an `unexpected-pass` case fail the run.

Exit codes: 0 when every selected case passed, 1 when a case failed or changed, 2 when the script could not complete (build failure, git error, invalid arguments, container failure).

### TreeTime v0 and Python scripts

TreeTime v0 (`packages/legacy/treetime`), the golden-master capture scripts, and the research scripts in `test_scripts/` share one Python environment, defined in `pixi.toml` and locked in `pixi.lock`. It pins every package and uses OpenBLAS 0.3.28 on a single thread, the same BLAS as the Rust build, so v0 and v1 run the same linear algebra kernels. The environment puts v0 of the checkout on `PYTHONPATH`, so both the v0 `treetime` command and `import treetime` in scripts work.

In the Python image:

```bash
./dev/docker/python treetime ancestral --help
./dev/docker/python python3 test_scripts/fitch.py
```

On the host, [pixi](https://pixi.sh) comes with `mise install` and creates the environment in `.pixi/` on first use:

```bash
pixi run treetime ancestral --help
pixi run python3 test_scripts/fitch.py
```

After changing `pixi.toml`, run `pixi lock` and commit both files. Package releases must be at least seven days old; `exclude-newer` in `pixi.toml` enforces this.

## Build profiles

- `dev` (`just build`, `just run`, the tests): unoptimized workspace crates, dependencies at `opt-level = 2`, full debug info
- `dev-opt` (`just run-dev-opt`): `dev` with optimized workspace crates, for long runs on real datasets; rebuilds are slower
- `release` (`just build-release`, `just run-release`, `just example`, `just smoke`): optimized and fast to rebuild, without LTO
- `dist` (`just build-dist`, `just cross`, `just build-napi dist`, `just cross-napi`, nightly releases): the shipped binary and Node addon, with fat LTO and one codegen unit
- `profiling` (`just build-profiling`, `just profile`): `dist` with full debug info
- `bench` (`just bench`): the `dist` settings

## Performance

- `just bench` runs the benchmarks
- `just build-dist treetime` builds the shipped binary into `.out/treetime`; use it with `hyperfine`, which the container provides
- `just profile treetime -- <args>` (host) samples a profile with samply or perf; read `dev/profile --help` first

## Dependencies

Every dependency release must be at least seven days old before the project adopts it:

- Bun refuses younger npm packages (`minimumReleaseAge` in `bunfig.toml`)
- Cargo has the age check as an unstable feature until Rust 1.100. `just deps-update` and `just deps-upgrade` enable it for their resolution, and `just deps-age` fails on a lockfile entry younger than seven days
- Python packages: `exclude-newer` in `pixi.toml` makes `pixi lock` skip younger releases
- Tools in `.config/mise.toml`: `minimum_release_age` makes `just tools-outdated` list only newer releases that are at least seven days old. mise does not filter an exact version, so check the date of a version you type by hand
- Base images in `dev/docker/` follow the same rule by hand
- mise itself: the images and CI install the version in `dev/docker/files/mise-version`, verified by its line in `dev/docker/files/checksums`. `min_version` in `.config/mise.toml` is the oldest mise that reads the configuration; raise it when the configuration needs a newer mise
- After changing a tool in `.config/mise.toml`, `just tools-lock <tool>` locks only that tool. Locking calls the GitHub API, which allows 60 anonymous requests per hour. With `MISE_GITHUB_TOKEN` set, mise authenticates instead, on the host and in the container alike (`dev/docker/run` forwards it), for example `MISE_GITHUB_TOKEN="$(gh auth token)" ./dev/docker/run just tools-lock <tool>`

Rust dependencies are pinned exactly in the workspace `Cargo.toml`, JavaScript dependencies exactly in the manifests, with shared packages in the Bun catalog of the root `package.json`. The React packages stay on 19.x, because Auspice runs in-process and each React major needs its own check of the Auspice modules TreeTime imports; `just lint-ts` enforces this.

The dependency recipes run in the main checkout only. `just audit` checks both dependency graphs against the security advisory databases.

## Cross-compilation and releases

`just cross` (host) builds the shipped CLI (`dist` profile) for every target in `dev/cross/targets` in its cross image, in parallel, into `.out/`. `just cross --target=aarch64-apple-darwin` builds one target. Without `just` on the host, run `./dev/cross/all treetime`.

`just cross-napi` (host) builds the Node addon of the desktop app the same way, with the same profile and CPU flags, for every target in `dev/cross/targets-desktop`, the targets Electron supports (no musl), into `.out/app-napi-<target>.node` (`./dev/cross/all --lib app-napi` without `just`). `dev/cross/check --lib` fails when an addon needs OpenBLAS, the Fortran runtime, a MinGW runtime DLL, or a macOS library outside the system. The Windows addon is built with MinGW, like the CLI, and loads into Electron, which is built with MSVC: Node-API reaches the addon through the functions that `napi-sys` looks up in the running executable, not through an import library.

The release binaries require an x86_64 CPU with AVX2 (Haswell or newer) and, on Linux aarch64, ARMv8.2.

### Nightly releases

Prerelease builds are published in [neherlab/treetime-nightly](https://github.com/neherlab/treetime-nightly/releases). A nightly publishes every CLI binary and desktop package that builds, without the linkage checks and start tests of the CI workflow.

1. `.github/workflows/schedule-nightly.yml` on `master` runs daily at 04:00 UTC and calls `.github/workflows/nightly.yml` on `rust`, which skips when `rust` has no new commits
2. `nightly.yml` builds the cross-compilation matrix of `.github/workflows/cli-build.yml`, and the desktop apps with `.github/workflows/desktop-build.yml`: the addon of each target in its cross image, the bundles, and the packages. A failed target does not block the others
3. `dev/publish-nightly` creates the prerelease with the built binaries and desktop packages. Its notes list the changelog since the previous nightly, followed by a downloads table that links every asset, one row per platform

`./dev/trigger-nightly` starts a nightly by hand. The version format is `<cargo-version>-nightly.<YYYYMMDD>T<HHMMSS>Z+<short-sha>`.

The same nightly also deploys the web app (next section). `./dev/trigger-nightly --only cli` runs only the CLI release, `--only desktop` only the desktop release, `--only web` only the web deploy; the scheduled run does all three. A desktop-only release has no CLI binaries, so the scheduled run still builds the CLI for the same commit, and a desktop-only run also runs when `rust` is unchanged.

### Nightly web deploy

The nightly deploys the web app (UI and API server) of `rust` to one Hetzner Cloud server, protected by one shared password. The server files are in `dev/deploy/hetzner/`.

How a deploy works:

1. The `deploy-web` job of `nightly.yml` calls `.github/workflows/web-deploy.yml` for the same commit as the CLI release. It runs when a CLI nightly is built, or on every `--only web` run, and does not depend on the CLI jobs
2. `dev/deploy/build-web-image <tag>` builds the static `treetime-server` binary for `x86_64-unknown-linux-musl`, the UI, and the image `treetime-web:<tag>` (`dev/deploy/web.dockerfile`), which only copies these artifacts and `data/`. It also builds `treetime-caddy:<tag>` (`dev/deploy/caddy.dockerfile`): Caddy with the rate limit plugin `github.com/mholt/caddy-ratelimit`, because stock Caddy cannot limit requests. The tag is the short commit hash. The same command builds the images locally
3. The workflow saves both images into one `image.tar` and sends it, `compose.yaml`, and `Caddyfile` as one zstd-compressed tar stream over SSH, without an image registry
4. On the server, the SSH key of user `deploy` runs only `/usr/local/bin/treetime-deploy` (a copy of `dev/deploy/hetzner/treetime-deploy`). It loads the image, installs the files in `/opt/treetime/`, sets `TREETIME_TAG` in `/opt/treetime/.env`, and restarts the Docker Compose stack. When the app does not answer within 60 s, it restores the previous release and the job fails. It keeps the images of the last three deployed releases
5. Caddy serves HTTPS with a Let's Encrypt certificate and asks for the password (HTTP Basic Auth) before it forwards requests to the app. Before the password check, it answers HTTP 429 to a client address (an IPv6 `/64` network counts as one address) above 600 requests per minute. Each new wrong password costs one bcrypt check, and the limit bounds the CPU that password guessing can take. The compose network has IPv6 enabled, because Docker otherwise forwards IPv6 connections through its userland proxy, and Caddy then sees the gateway address of the network in place of the client address. Caddy serves HTTP/1.1 and HTTP/2 without HTTP/3, because QUIC costs more CPU than TCP on the two cores the runs need, sends `Strict-Transport-Security`, and leaves compression to the app. The app answers only the host `app` of the deploy health check and `TREETIME_DOMAIN` (`--allowed-host` in `compose.yaml`)

The server has no login of its own, lists all runs to everyone who has the password, has no limit on parallel runs, and never deletes runs: they stay in `/opt/treetime/runs/`. The container limits of `compose.yaml` (2 CPUs, 3 GB memory, 50 MB of uploads per run) bound what one client can take from the server.

A failed scheduled deploy is retried only when `rust` changes. Run it again with `./dev/trigger-nightly --only web`.

#### One-time setup

1. **Server**: create a Hetzner Cloud CX23 (x86, Ubuntu 24.04) with the admin SSH key and the user data `dev/deploy/hetzner/cloud-init.yaml`, after replacing `<deploy public key>` in it (step 2). Attach a Hetzner Cloud Firewall that allows inbound TCP 22, 80, and 443 only. Check that `grep -c avx2 /proc/cpuinfo` prints a non-zero count: the binary requires AVX2
2. **Deploy access**: create a key pair for the deploy, for example `ssh-keygen -t ed25519 -N '' -C treetime-deploy -f deploy_key`, and put the public key in the user data (step 1), or later in `/home/deploy/.ssh/authorized_keys` with the same `restrict,command=...` prefix. As admin, copy `dev/deploy/hetzner/treetime-deploy` to `/usr/local/bin/treetime-deploy` with mode 755; copy it the same way after each change of the script. In the GitHub repository, create the environment `nightly-web`, limit its deployment branches to `rust` and `master`, and set there:
   - secret `DEPLOY_SSH_KEY`: the private key
   - secret `DEPLOY_KNOWN_HOSTS`: the output of `ssh-keyscan <host>`
   - variable `DEPLOY_HOST`: the host name
3. **DNS**: before the first deploy, add `A` and `AAAA` records in Route 53 for the host name and set `TREETIME_DOMAIN=<host name>` in `/opt/treetime/.env`. Caddy requests the certificate when it starts, and Let's Encrypt rate-limits repeated failures, so the records must resolve first
4. **Password**: before the first deploy, create the hash with `docker run --rm -it caddy:2 caddy hash-password --bcrypt-cost 8` and set `TREETIME_USER=<user name>` and `TREETIME_PASSWORD_HASH='<hash>'` in `/opt/treetime/.env`. Keep the single quotes: the hash contains `$`, which Docker Compose expands in unquoted and double-quoted values. Setting a new hash revokes all access; share the new password with the people who need it. The password only keeps crawlers and bots away, so the hash uses bcrypt cost 8 in place of Caddy's default 14: on a CX23, a wrong password then costs about 20 ms of CPU in place of about 1.4 s

`/opt/treetime/.env` belongs to user `deploy` with mode 600. The deploy changes only its `TREETIME_TAG` line.

#### Trust model

Membership in the `docker` group equals root on the server, and each deploy sends its own `compose.yaml`, which can mount any host path. Whoever holds the deploy key therefore controls the server. This is accepted because the server holds nothing secret beyond its TLS certificate and the password hash, and only workflows on `rust` and `master` can read the key (environment `nightly-web`).

#### Rollback and maintenance

To run an older release, as admin set `TREETIME_TAG` in `/opt/treetime/.env` to one of the tags of `docker image ls treetime-web` (the image `treetime-caddy` has the same tags), then run `docker compose up -d` in `/opt/treetime/`. The next deploy replaces it.

The server installs security updates by itself and reboots at 03:00 UTC when an update requires it. Docker starts on boot and brings the stack back.

## Continuous integration

`.github/workflows/cli.yml` runs on pull requests and on pushes to `rust`:

- The check groups of `just check-all` in parallel jobs, except `hawk`, which needs a second full compilation; run it locally
- On pushes to `rust`: the release builds of every target and compatibility runs on Linux distributions, macOS, and Windows

The CI jobs pull the container images from Docker Hub by the hash of their build inputs, and build them when the inputs changed. Only pushes to `rust` publish images.
