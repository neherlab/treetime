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
- The images link OpenBLAS, liblzma, and libzstd as static prebuilt archives instead of compiling them in every new target directory: OpenBLAS of [binarylandia/build_openblas](https://github.com/binarylandia/build_openblas), built from one OpenBLAS develop commit with clang and flang of LLVM 22, with the Fortran runtime of flang (`libflang_rt.runtime.a`) that its LAPACK calls, and liblzma and libzstd of [binarylandia/staticbourg](https://github.com/binarylandia/staticbourg). `.config/mise.image.toml` pins one archive of each library per target, and `.config/mise.image.lock` locks them by URL and checksum. `dev/docker/files/install-image-libraries` installs the archives of the image's targets, the Windows targets included, and writes their pkg-config folders into the `[env]` table of `/.cargo/config.toml`, so the build scripts of `treetime-openblas`, `liblzma-sys`, and `zstd-sys` link the archives through pkg-config; the threaded WebAssembly build takes its own OpenBLAS from `/.cargo/wasm-mt-env.toml`. The folders are named after the archive versions, so a new archive changes the environment of the build scripts, and cargo reruns them and relinks every binary, test, and addon
- The image build stops when the build record of an archive names another LLVM release, target toolchain, or Emscripten version than the image's entries, or other CPU flags than `export_dist_flags` of `dev/lib/dist-flags.sh`. `just prebuilt-libs-check` fails when the liblzma and libzstd archives have other C versions than `liblzma-sys` and `zstd-sys` of `Cargo.lock` bundle (a crate update then needs a new staticbourg release first), when the entries of one library name different releases, and when the CPU level of the development builds in `.cargo/config.toml` differs from that of `export_dist_flags`
- The native image compiles and links x86_64 code with clang of LLVM 22 on the sysroot, crt files, libgcc, and libstdc++ of the GCC 16.2.0 toolchain of the shipped Linux build, so the binaries built there need glibc 2.17 or later and also run on the host. Its `/.cargo/config.toml` adds mold as the linker of development builds
- Commands run as your user, with all Linux capabilities dropped
- Each container is named `treetime-<image>-<start time>-<process id>`, such as `treetime-native-20261022-011500-4242` or `treetime-aarch64-apple-darwin-20261022-011500-4250`, so `docker ps` shows what it runs. It uses at most `TREETIME_DOCKER_CPUS` CPUs (see [Machine settings](#machine-settings))
- On Apple Silicon and other arm64 hosts the image runs under x86_64 emulation, which is slower
- The container uses host networking, so the dev servers are reachable from the host browser. On Docker Desktop, turn on host networking in the settings

### Host

Requirements: [mise](https://mise.jdx.dev), [rustup](https://rustup.rs), bash 4 or later (on macOS: `brew install bash`), and a C toolchain with pkg-config. On Debian and Ubuntu:

```bash
sudo apt-get install build-essential pkg-config
mise install
just setup
just check
```

Host builds go to `.build/host/`. Host builds of the Rust code are tested on Linux x86_64 only; on macOS, build in the container. The dylint and IQ-TREE tools are available on Linux only; run the recipes that need them in the container on other hosts.

With mise active, `mise install` also installs the static OpenBLAS and liblzma archives of the shipped Linux and macOS builds (the `openblas` and `liblzma` entries of `.config/mise.toml`), and `HOST_PKG_CONFIG_PATH` puts their pkg-config folders before your own `PKG_CONFIG_PATH`, so native builds link these archives statically, as the release does. libzstd compiles from the bundled C sources of `zstd-sys`.

Without the mise archives, a build links the OpenBLAS of the system, which the build script of `treetime-openblas` finds through pkg-config; when it finds none, the build stops with the install commands of each platform (`packages/treetime-openblas/build.rs`). On Debian and Ubuntu, install `libopenblas-dev`; pkg-config finds it without further settings. Such a build also needs `protoc` for the UShER schemas (`protobuf-compiler` on Debian and Ubuntu), which mise otherwise provides. A system OpenBLAS differs from the project builds:

- A multi-threaded OpenBLAS, such as the default `libopenblas-pthread-dev` of Debian and Ubuntu, starts its own threads next to the threads of TreeTime, which compete for the CPUs. `OPENBLAS_NUM_THREADS=1` avoids this
- Its results can differ in the last digits from those of the project builds

Remove any `[target.<target>.openblas]` table from `~/.cargo/config.toml` that an earlier host setup added: such an entry replaces the build script of `treetime-openblas` without a message.

Host builds link with rust-lld, the default linker of rustc. mold is optional: add `"-Clink-arg=-fuse-ld=mold"` to the `rustflags` of `[target.x86_64-unknown-linux-gnu]` in `~/.cargo/config.toml`. A relink of the CLI after a one-file edit takes about as long with mold as with rust-lld, because the compile of the changed crate dominates.

### Machine settings

Optional settings go into the gitignored `.env` in the checkout; `.env.example` lists them:

- `KACHE_STORE`: directory of the [kache](https://github.com/kunobi-ninja/kache) compiler cache stores. Every build and clippy pass then compiles through kache, which shares compiled crates across the worktrees of this project. Workspace crates of incremental builds (`dev`, tests, `release`, clippy) keep their incremental state, so an edit rebuilds as fast as without kache; kache serves their dependencies and every build without incremental state (CI, `dist`, `profiling`, `bench`, cross builds). Dylint, hawk, and coverage compile without kache. The directory holds one store per kache version, environment (`host` or `docker`), and pass (`build`, `clippy`, or `cross-<target>`)
- `KACHE_MAX_SIZE`: size limit of each kache store, 100 GiB when unset
- `TREETIME_DOCKER_CPUS`: CPU budget of the containers that the dev scripts start, 8 when unset and at most the CPUs of the Docker daemon. A container of `dev/docker/run` uses at most this many CPUs; `dev/cross/all` and the Linux distribution tests divide the budget between the containers they run at the same time. Cargo in a container runs one job per whole CPU of its budget unless `CARGO_BUILD_JOBS` is set, so the compiler work stays within the budget also where a host daemon such as ananicy moves compilers out of the container cgroup
- `SMOKE_STORE`: absolute directory of a smoke snapshot store shared by all checkouts of this project on the machine, so each baseline is built once (see [Snapshot store](#snapshot-store)). Unset: each checkout keeps its own `snapshots/`
- `SMOKE_MIN_FREE_GIB`: free disk space that `dev/smoke` keeps on the file systems of the checkout and the snapshot store, 20 GiB when unset (see [Disk limits](#disk-limits))

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
| Clippy of the WebAssembly back end               | `just lint-wasm`, `just lint-wasm-mt` |
| Smoke test of the WebAssembly back end, in Bun   | `just test-wasm`                     |
| Smoke test of the threaded WebAssembly back end  | `just test-wasm-mt`                  |
| TypeScript lints                                 | `just lint-ts`                       |
| Format one toolchain                             | `just fmt-rs`, `fmt-ts`, `fmt-other` |
| Lints and tests of the tools in `dev/`           | `just check-dev`                     |

In the container, prefix each command with `./dev/docker/run`.

Recipe names follow one scheme:

- **Leaf recipes** run one tool in one mode. The suffix `-rs`, `-wasm`, or `-ts` names the toolchain: Rust, Rust for WebAssembly, or TypeScript, and `-wasm-mt` the threaded WebAssembly build; tools that serve one toolchain only keep their own name, such as `dylint`, `hawk`, or `knip`
- **Combined recipes** without a language suffix, such as `lint`, `test`, `fmt`, or `fix`, run the leaves of both toolchains and call no tool themselves
- **The suffix `-all`** adds the slow tools to the fast set of the same name: `fix-all` adds the dylint fixes, `lint-all` the lint libraries, hawk, and the dependency and config lints, and `check-all` the full gate

While working on one toolchain, run its leaves (`just fmt-rs`, `just l`, `just t`) and leave the slow tools to `check-all`.

The lint libraries and oxlint plugins in `dev/lints` have their own recipes, which `check` and `check-all` leave out: `just check-dev` runs `lint-dev-rs` and `test-dev-rs` for the dylint libraries, each in its own workspace with its pinned nightly, `lint-dev-ts` and `test-dev-ts` for the oxlint plugins, and `test-dev-sh` for the shell libraries in `dev/lib` (the port slots of `dev/lib/ports.sh`). Run them after a change in `dev/lints` or `dev/lib`.

`check` and `check-all` run their checks in parallel through `dev/run-checks`, keep going past failures, and list the failed checks at the end. Each check writes its output to `tmp/checks/<check>.log`. Warnings fail the checks, and a missing tool is a failure, not a skipped check. The full gate consists of groups, one CI job each; `just check-group <group>` runs one group serially, as its CI job does.

### Examples

Replace `$v` with a dataset path under `data/`, for example `flu/h3n2/20`:

```bash
just run treetime ancestral --method-anc=marginal --tree=data/$v/tree.nwk --alignment=data/$v/aln.fasta.xz --output-all=tmp/ancestral/$v
just run treetime clock --tree=data/$v/tree.nwk --metadata=data/$v/metadata.tsv --output-all=tmp/clock/$v
just run treetime timetree --tree=data/$v/tree.nwk --metadata=data/$v/metadata.tsv --alignment=data/$v/aln.fasta.xz --output-all=tmp/timetree/$v
```

## Lint settings

The custom dylint library and the oxlint plugin of `dev/lints/oxlint/` are both named `custom`, so a suppression reads `#[cfg_attr(dylint_lib = "custom", expect(no_comments, reason = ".."))]` or `// oxlint-disable-next-line custom/no-vague-identifiers -- reason`.

Clippy rejects `#[allow]` on an item: suppress a lint on the narrowest item with `#[expect(<lint>, reason = "..")]`, which fails once the lint no longer fires, so a stale suppression cannot stay.

The project settings of the lints are in these files:

- `justfile` variables: `dylint_rustflags` (the lint levels of the custom library), the check groups `checks_*` with `check_fast` and `check_full`, `lint_fast`, `lint_full`, `public_api_crates`, `dockerfiles`, and `react_pin_reason`
- `dylint.toml`: the settings of the dylint libraries, such as the render sources of `no_comments`, the entry points of `forbidden-reach`, the error style of `proper_error_type`, and the helper macros of `prefer_error_macros`
- `oxlint.config.ts`: a call of `projectConfig` from `dev/lints/oxlint/config.ts` with the package layout, the web scopes, the import and property restrictions, and the contracts package
- `clippy.toml` (thresholds and banned types and functions; `dev/lints/clippy.toml` replaces it for the lint libraries, with the same thresholds and without the bans), `oxfmt.config.ts`, `taplo.toml`, `.config/knip.json`, `.config/deny.toml`, `.config/hawk.toml`, `.config/jscpd.json`, and `.config/mordant-baseline.toml`

Rust comments are banned in TreeTime: every lint of the custom library is on, so `no_comments` keeps only the doc comments that a render source of `dylint.toml` carries into user-facing output, and comments that start with a marker such as `SAFETY:` or `TODO`. TypeScript comments are banned by the oxlint base configuration `dev/lints/oxlint/config.ts`, apart from tool directives.

## Project structure

```
packages/
  app-backend/       Run records, results, settings, and pre-run checks of the server, the addon, and the web2 worker
  app-cli/           CLI binary (treetime)
  app-commands/      Command configs and runners shared by the CLI, the server, and the addon
  app-contracts/     TypeScript types and schemas of the bridge between the apps and the Rust code
  app-datasets/      Bundled example datasets
  app-desktop/       Electron shell of the desktop app
  app-docs/          Docs site (Astro with headless fumadocs-core)
  app-napi/          Node addon of the desktop app (Rust, napi-rs)
  app-output/        Output writers shared by the CLI, the server, and the addon
  app-server/        HTTP API server of the web app (treetime-server)
  app-ui/            React components shared by the desktop and web apps
  app-wasm/          WebAssembly back end of the web2 app (Rust, Emscripten, wasm-bindgen)
  app-web/           Web app (Vite)
  app-web2/          Web app on the WebAssembly back end in a web worker (Vite)
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

The TypeScript apps (`app-ui`, `app-web`, `app-web2`, `app-desktop`, and the docs site `app-docs`) are a user interface around the Rust core. Algorithms, file parsers, readers and writers, configuration handling, and domain rules live in Rust and are tested there. The apps render what the Rust operations return and keep only presentation logic: layout, interaction state, and formatting for display. One implementation of each rule keeps the CLI, the web app, and the desktop app from diverging, and one test suite covers it.

- **CLI knowledge**: flags, value syntax, setting groups, defaults, the equivalent command line, and the YAML config come from clap and the config types. The setting catalog (`packages/app-contracts/setting-catalog.json`, generated by `just gen catalog` and imported from `@neherlab/app-contracts/catalog`) describes every command and setting: titles, main settings, and the settings that act as one choice. `check-config` reports which option of each choice a configuration selects, and `run-config` and `check-config` render the command line and YAML
- **Checks**: `check-config` classifies configuration problems and input facts into checks that block the run, warn, or advise, with the settings that fix them
- **Types**: every value that crosses from Rust to TypeScript is a Rust type with a derived schema in the OpenAPI document; the apps read it through the generated types and zod schemas in `app-contracts`, never through hand-written ones

### Interface components

The interface is built from [shadcn/ui](https://ui.shadcn.com) components in the `base-vega` style, which run on Base UI. The components are vendored into `packages/app-ui/src/ui/`, one file per component, and the app composes them; it does not style raw elements. TreeTime owns the vendored files: they pass the same lint rules as the rest of the source and keep only the parts the app uses. `packages/app-ui/src/ui/UPSTREAM.md` records their source, the rewrites that turn a registry file into a component, and the local changes to re-apply when a component is refreshed.

The theme tokens (colors, radius, fonts) live in `packages/app-ui/src/tokens.css`: green-grey surfaces with a deep teal accent, Lato for text, IBM Plex Sans for headings, and IBM Plex Mono for code. `packages/app-ui/src/components.css` holds the imports that the components need, `packages/app-ui/src/base.css` the base rules that every surface shares, and `packages/app-ui/src/theme.css`, the stylesheet of the app, imports the three and adds the page rules of the app. `packages/app-ui/src/ui/shadcn.css` is the stylesheet of the `shadcn` npm package that defines the variants the components use.

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

Each checkout takes its ports from its port slot, a group of 10 ports in the block 14000-14999, so worktrees run side by side. The main checkout holds slot 0 (ports 14000-14007). A linked worktree claims the lowest slot that no other checkout of the repository holds and whose ports have no listener, and keeps it until the worktree is removed; the claims are files in the git common directory (`dev/lib/ports.sh`). `just status` prints the slot and its ports. `TREETIME_API_PORT`, `TREETIME_WEB_PORT`, `TREETIME_SERVE_PORT`, and the other `TREETIME_*_PORT` variables override single ports.

A route launcher on the host can publish the web dev server through a local HTTPS proxy, such as [portless](https://github.com/vercel-labs/portless), under a hostname such as `https://treetime.localhost`. The launcher starts `just up` as its child process, assigns the hostnames, and removes the route when the command ends. The checkout and the launcher exchange two files in `tmp/portless/`:

- `web.hostnames.json`: `{"hostnames": [...], "urls": [...]}`, written by the launcher before it starts the command. `just status` prints the URLs. Without this file, the server has no public hostname
- `web.json`: `{"port": <n>}`, written by `just up` after the web server answers, through a temporary file and a rename, so the launcher never reads a partial file. `just up` writes it also when no launcher runs

The checkout never reads or writes the state of the proxy itself. The Vite server accepts every `*.localhost` hostname, so each assigned name works without configuration.

#### Storage

The web app has no database. The server keeps each run in its own folder `<runs-dir>/<run-id>/`: the run record `run.json`, the event log `events.jsonl`, the uploaded files in `inputs/`, and the results in `out/`. The browser keeps the sidebar width in `localStorage` under the key `treetime-preferences`, so every tab shares it, and the unfinished form in `sessionStorage` under the key `treetime-draft`, so each tab keeps its own form through a reload. The theme stays with `next-themes`, which stores it in `localStorage` under the key `theme` and applies a change in one tab to the other tabs at once. The server reads run inputs only from the `inputs/` folders of the runs; the web app uploads a chosen example file into the run like a user file (see [Example data](#example-data)).

Each app mode keeps its runs in its own directory:

- `tmp/app/web-dev/runs`: development servers (`just up`)
- `tmp/app/web-prod/runs`: production server (`just serve`)
- `tmp/app/treetime-dev/runs`: desktop app in development mode (`just desktop`); see [Desktop app](#desktop-app) for its other folders
- `tmp/app/treetime-prod/runs`: production build of the desktop app started from the checkout (`just desktop-prod`)
- `<app folder>/runs`: installed desktop app, unless the settings or `TREETIME_RUNS_DIR` name another folder; see [Desktop app](#desktop-app)

A deployment passes `--runs-dir` to `treetime-server` itself.

#### Event streams and HTTP/2

The web app follows changes through server-sent event streams: `GET /api/events` for changes to any run, and `GET /api/runs/{id}/events` for each run it follows. Every open stream holds one connection. Over HTTP/1.1, browsers allow at most 6 connections per host, shared by all tabs, so more open streams stall every further request to the server. Over HTTP/2, all requests share one connection.

The API server speaks HTTP/1.1 without TLS, and browsers use HTTP/2 only over TLS. A deployment therefore puts a reverse proxy in front of the server that terminates TLS and serves HTTP/2 to the browser, for example nginx with `listen 443 ssl; http2 on;`. In development, the HTTPS proxy of a route launcher (see [Web app](#web-app)) does this. Without such a proxy, at most 6 streams stay open per host across all tabs.

#### Caching, compression, and allowed hosts

`treetime-server` serves the UI build (`STATIC_DIR`) and the API with these rules, so a proxy in front of it needs no caching or compression settings of its own:

- **UI files**: the file names in `assets/` contain a hash of their content, so they are cached for a year (`Cache-Control: public, max-age=31536000, immutable`). Every other file, `index.html` included, is revalidated on each use (`no-cache`), so a new deploy reaches the browser at the next page load. Only page navigations (`Accept: text/html`) fall back to `index.html`; a missing asset or API path answers 404
- **Compression of UI files**: `bun run build:web` writes a gzip (level 9) and a brotli (quality 9) copy next to each text file of the build, and the server sends the copy the browser accepts, with `Vary: Accept-Encoding`. Brotli quality 11 would make the files about 8 % smaller but adds seconds to every build
- **Compression of API responses**: zstd, or gzip for browsers without zstd, at the fastest level. Run outputs of large alignments shrink far more with zstd than with gzip, because gzip only finds repeats within 32 KB, less than one genome. Event streams and already compressed files (`.xz`, `.gz`, `.zst`, fonts) stay uncompressed
- **API caching**: API responses are not stored (`no-store`). The results, the Auspice document, and the comparison of finished runs carry an `ETag`, so a browser that asks again gets `304 Not Modified` without the server rebuilding the answer. The tag changes when the server restarts
- **Downloads**: output files and run archives stream to the browser, so the server holds only a small buffer per download. Archives store their files uncompressed and rely on the HTTP compression above
- **Allowed hosts**: the server answers only requests for `localhost`, `*.localhost`, `127.0.0.1`, `[::1]`, and the hosts given with `--allowed-host` (repeatable) or `ALLOWED_HOSTS` (comma-separated); other hosts get HTTP 403. This stops other websites from reaching a local server through DNS rebinding. The server sends no CORS headers, because the UI calls it from its own origin

### Web2 app

The web2 app (`packages/app-web2`) is the web app without a server: the interface of `app-ui` sends its API requests to the WebAssembly build of the API router (`packages/app-wasm`), which runs in a dedicated web worker of the page. [kb/decisions/app-web2-wasm-worker-back-end.md](../../kb/decisions/app-web2-wasm-worker-back-end.md) describes the design.

| Task                                                                       | Command                                   |
| -------------------------------------------------------------------------- | ----------------------------------------- |
| Build the single-threaded WebAssembly module into `packages/app-wasm/pkg/st/` | `just build-wasm [dev\|release\|prod]`    |
| Build the threaded WebAssembly module into `packages/app-wasm/pkg/mt/`     | `just build-wasm-mt [dev\|release\|prod]` |
| Clippy of the single-threaded and of the threaded module                   | `just lint-wasm`, `just lint-wasm-mt`     |
| Smoke test of the single-threaded module, in Bun                           | `just test-wasm`                          |
| Smoke test of the threaded module with 2 and 8 threads, in Node            | `just test-wasm-mt`                       |
| Start the web2 dev server                                                  | `just up-web2`                            |
| Build web2 for production and serve it                                     | `just serve-web2`                         |
| Production build of web2 into `packages/app-web2/dist/`                    | `just build-web2`                         |

`just build-wasm` builds `app-wasm` for `wasm32-unknown-emscripten` with cargo. rustc links the module through `emcc` of Emscripten, and `emcc` runs the `wasm-bindgen` CLI on the linked module, so the build writes one JavaScript module with the Emscripten runtime and the bindings, the WebAssembly binary, and the TypeScript declarations; the recipe copies them into the package `@neherlab/app-wasm` (`packages/app-wasm/pkg/`). Emscripten provides the C library of OpenBLAS and of the compression libraries, an in-memory file system for the runs, clocks, and entropy. The dev image installs Emscripten, a single-threaded and a threaded OpenBLAS that [binarylandia/build_openblas](https://github.com/binarylandia/build_openblas) builds with `emcc` for this Emscripten version, and the WebAssembly liblzma and libzstd of staticbourg, all from `.config/mise.image.toml` (the `http:emscripten` and `*-wasm32-unknown-emscripten` entries). Emscripten gives no ABI guarantee between versions, so the image build fails when an archive comes from another Emscripten version, and an Emscripten upgrade needs new `build_openblas` and staticbourg releases first. The `prod` mode uses the `dist` profile, which `emcc` optimizes at link time. `just gen wasm-types` writes the declarations `packages/app-wasm/pkg/app_wasm.d.ts`, the only committed file of the package, so the TypeScript checks run without a WebAssembly build. The version of the `wasm-bindgen` CLI in `.config/mise.toml` and the version of the `wasm-bindgen` crate in `Cargo.toml` must be equal, because the CLI accepts only modules of its own schema version. mise installs the CLI from [binarylandia/appville](https://github.com/binarylandia/appville) (task `tasks/wasm-bindgen`): upstream wasm-bindgen with a patch that keeps the start function of the linker, which the threaded module needs to set up its shared memory. A new crate version needs a new appville release first.

The common link settings of both modules are in `[target.wasm32-unknown-emscripten]` of `.cargo/config.toml`, among them `-sSTACK_OVERFLOW_CHECK=2` and `-Wl,--stack-first`, which stop the module with "stack overflow" instead of letting an overflow overwrite data. Each cargo call for the target adds the settings of its build with `--config`: `dev/cargo/wasm-st.toml` for the single-threaded module and `dev/cargo/wasm-mt.toml` for the threaded one. `just build-wasm-mt` builds on the nightly of `dev/docker/files/wasm-threads-toolchain`, which rebuilds the standard library with atomics (`-Zbuild-std`), in its own target folder, and links with `-pthread` and the threaded OpenBLAS. The module starts a pool of `jobs + 2` pthread workers with the module, because rayon waits for every thread of its pool and a new pthread worker starts only after its creator returns to the event loop. `just build-wasm-mt` fails when the threaded module declares other exports than `pkg/app_wasm.d.ts`. The package exports the single-threaded module as `@neherlab/app-wasm` and the threaded one as `@neherlab/app-wasm/threads`.

`just lint-wasm` and `just lint-wasm-mt` run clippy over the binary of `app-wasm` for the target, the second over the code that compiles only with atomics; `just check` runs both. The tests of `app-wasm` run natively, because a test harness needs unwinding and the module aborts on a panic. `just test-wasm` builds the single-threaded module and runs `packages/app-wasm/smoke.ts` in Bun: it creates a back-end instance and a job instance of the module and drives the job protocol without workers, with an `ancestral` run of `data/flu/h3n2/20` and its files and archive, a request answered while a job computes, a cancelled run, and an input that fails its checksum. `just test-wasm-mt` runs the same smoke test on the threaded module with 2 and with 8 threads under Node, because Bun does not start the pthread workers of Emscripten, and Node stays alive during a threaded run only through a timer that the test holds; it is slow and outside every gate. `just up-web2` builds both modules once; run `just build-wasm` and `just build-wasm-mt` again after a Rust change. `TREETIME_WEB2_PORT` and `TREETIME_WEB2_SERVE_PORT` override the ports of the dev server and of the served production build (`just status` prints them).

The page compiles the single-threaded module once and starts the back-end worker with it (`WorkerConnection`, `packages/app-web2/src/backend/connection.ts`). The worker answers each request through an exchange of the module and sends the response head and body chunks back through a message port (`ExchangeServer`, `packages/app-web2/src/backend/serve.ts`). The page uses the fetch of the desktop renderer for this port (`createPortFetch()` of `@neherlab/app-ui/port-fetch`), and transfers upload bodies instead of copying them. A back-end worker that fails to start stops the back end; a back-end worker that crashes restarts with the same module, and the requests in flight and every run of the page are lost. The content security policy of the web2 page allows `'wasm-unsafe-eval'`, which WebAssembly compilation needs, and nothing else that evaluates code.

Each run computes in a job worker of its own, which the back-end worker starts for the run (`JobHost`, `packages/app-web2/src/backend/jobs.ts`, and `packages/app-web2/src/job/worker.ts`), so several runs can compute at the same time. The back-end worker sends the job its input files, writes each progress event of the job into the run's event log, so the page shows the progress live, and receives the output files. Stop terminates the job worker at once. A panic, an out-of-memory abort, or a stack overflow in a job fails only that run, with the stderr lines and the abort reason of the module as the error. The job worker runs the threaded module, 8 threads by default; the URL parameter `jobs` sets the number of threads of each run (`?jobs=4`), and `?jobs=1` selects the single-threaded module. Threads need a cross-origin isolated page in a secure context (HTTPS or `localhost`): the dev server and the preview server of `just serve-web2` send `Cross-Origin-Opener-Policy: same-origin` and `Cross-Origin-Embedder-Policy: require-corp` (`packages/app-web2/vite.config.ts`). Without isolation, without `SharedArrayBuffer`, or when the browser cannot create a shared memory of 4 GiB, and after a threaded job fails to start, the jobs run the single-threaded module, and the console says why.

The back end has the routes of `treetime-server`, without the routes of a local app (settings, workspace, examples download). The runs, the uploads, and the outputs live in the in-memory file system of the back-end worker and stay in the tab until it reloads or closes; the inputs section of the page says so, and web2 downloads run files and archives through the generated client. The request timeout and the keep-alive of the event streams exist only in the native builds. [kb/decisions/app-web2-wasm-worker-back-end.md](../../kb/decisions/app-web2-wasm-worker-back-end.md) lists the open issues.

### Docs site

The docs site (`packages/app-docs`) is a prerendered multi-page site in the look of the app. It is its own [Astro](https://astro.build) site: the headless library [`fumadocs-core`](https://www.fumadocs.dev) supplies the page tree, the language fallback, the neighbour links, and the tracking of the visible heading, and `@neherlab/app-ui` renders every visible part. The docs are a separate build, not a route of the app, because they change once per release while the app deploys on every commit.

| Task                                                                       | Command           |
| -------------------------------------------------------------------------- | ----------------- |
| Start the dev server for writing pages, with live reload                   | `just docs`       |
| Build the docs as published and serve the build like the docs host         | `just docs-prod`  |
| Build into `packages/app-docs/dist/` and `packages/app-docs/dist-root/`    | `just build-docs` |

`just docs` runs in the `sess` session `tt-dev-docs` and `just docs-prod` in `tt-dev-docs-prod`. The dev server serves `http://127.0.0.1:<port>/dev` on the docs port of the checkout (14005 in the main checkout) and reloads pages and styles as they change; it runs no post-build step, so search and the version list are missing there. `just docs-prod` builds the site and serves it on the docs serve port (14006 in the main checkout) with the file layout of the docs host: the site root files (`versions.json`, the page that leads to the preferred version) at `/`, the build at `/dev`, pages without their `.html` extension, and the not-found page of the locale for a missing page (`packages/app-docs/src/preview/serve-site.ts`). `TREETIME_DOCS_PORT` and `TREETIME_DOCS_SERVE_PORT` override the ports, and `just status` prints them.

The app and the docs share three layers, and each surface wires the rest itself:

- **Styles**: `packages/app-docs/src/styles/docs.css` loads Tailwind, `tokens.css`, `components.css`, and `base.css` of app-ui, and adds the prose styles of `@tailwindcss/typography` with the colors of the app tokens. An app-ui component inside an article sits in a `not-prose` element, so it looks as it does in the app
- **Components**: the header parts (`packages/app-ui/src/chrome/`), the theme control, the search palette, the sidebar items, the alerts, and the setting presentation (`packages/app-ui/src/settings/help.tsx`) are app-ui components. The lint bans `@base-ui/react`, `react-aria-components`, and `next-themes` in the docs, so a missing control goes into `packages/app-ui/src/ui/`
- **Messages**: the strings of the shared parts and of the docs layout are Paraglide messages of app-ui (see UI messages below)
- **Surface wiring**: the app keeps its router, API, and stores; the docs wire the content collections, `fumadocs-core`, Pagefind, and the version data

The layout (`src/layouts/DocsLayout.astro`) has the header (`src/components/DocsHeader.astro`), the navigation tree, the article, and the table of contents. Astro renders each React component as its own tree, and a `client:` directive works only on a component that the `.astro` or `.mdx` file imports itself, so every interactive part is one island: the header actions (search, version select, theme), the mobile navigation, the table of contents, the version banner, and the interactive author components. Everything else renders at build time and ships no JavaScript. `src/lib/urls.ts` builds every link of the site. Logic lives in `.ts` and `.tsx` files, and `.astro` files stay thin, because the type check does not read `.astro` files.

Content conventions:

- Pages are `.mdx` files in `src/content/docs/`, German translations in `src/content/docs/de/`; `meta.json` files order the navigation. A German path without a German file shows the English page with a notice
- Callouts are directives (`:::note`, `:::tip`, `:::warning`, `:::danger`); another type stops the build. A heading that starts with `1. ` becomes a step
- Links are root-relative without the version (`/showcase/math`); the Markdown pipeline (`src/markdown/pipeline.ts`) adds the version base and the language prefix
- Math uses KaTeX with `\tag{n}` and an HTML anchor before the display block, because KaTeX has no `\label`; shared macros are in `src/math/macros.ts`
- Interactive components (`SyncedTabs`, `Disclosure`) need an import and a `client:` directive in the page; `FileTree`, `LinkCard`, and `CardGrid` need neither
- The settings pages (`src/content/docs/reference/settings/`) render the setting catalog, one page per command, and the build fails when a command has no page; the command-line reference (`reference/cli.md`) is generated by `just gen cli-docs`. Never edit generated text by hand
- The showcase pages hold placeholder text that exercises each feature; real guides replace them. `showcase/writing.mdx` lists the syntax of each feature

`astro build` alone produces the complete output: one post-build step (`src/integrations/post-build.ts`) writes the Pagefind search index under the canonical page URLs, the page manifest `pages.json`, the site root files in `dist-root/`, and the `.gz` copies, and checks that each executable inline script of each page has its hash in the content security policy.

#### UI messages

The UI strings of the shared parts, the docs, and (in English) the app are Paraglide messages: `packages/app-ui/messages/en.json` and `de.json`, in the inlang project `packages/app-ui/project.inlang`, whose `settings.json` loads its plugins from `node_modules`. `just gen i18n` compiles them into `packages/app-ui/src/i18n/generated/`, which the builds, the tests, and the type check read; no build runs the Paraglide Vite plugin. To add a message, add the key to both files and run `just gen i18n`; a test fails when a key or a parameter is missing in a locale. A component reads its locale with `useLocale()` and passes it to each call, `m.<key>(params, { locale })`, so pages that render side by side at build time cannot mix languages. A new locale goes into `settings.json`, a message file, and the `i18n` and route settings of the docs. Setting texts, help texts, CLI help, errors, and run warnings come from Rust and stay English on every page.

### Desktop app

`just desktop` starts the Electron app in development mode, with the Vite server of its renderer on the desktop port of the checkout (14007 in the main checkout, `TREETIME_DESKTOP_PORT` overrides it). `just desktop-prod` builds the desktop app for production (the Node addon in the `dist` profile and the bundled renderer) and starts that build from the checkout. In the container both need the host display: `TREETIME_DOCKER_X11=1 ./dev/docker/run just desktop`.

The Node addon (`packages/app-napi`, a Rust `cdylib`) is built by cargo, like the CLI, and copied to `packages/app-napi/app-napi.node`, the `main` file of the package: `just build-napi` builds the `dev` profile for `just desktop`, and `just build-napi dist` the `dist` profile for `just build-desktop` and `just desktop-prod`. `just gen napi-types` generates its TypeScript types (`index.d.ts`).

#### Folders and settings

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

#### Example data

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

#### Packages

`just package-desktop [target]` packages the desktop app with [electron-builder](https://www.electron.build) into `.out/treetime-desktop-<target>.<ext>`:

- Linux (`x86_64-unknown-linux-gnu`, `aarch64-unknown-linux-gnu`): an AppImage, which runs without installation
- Windows (`x86_64-pc-windows-gnullvm`, `aarch64-pc-windows-gnullvm`): a zip archive; users extract it and start `TreeTime.exe`
- macOS (`x86_64-apple-darwin`, `aarch64-apple-darwin`): a dmg

The recipe builds the addon for the native target `x86_64-unknown-linux-gnu` itself; the other targets take the addon from `just cross-napi` (see [Cross-compilation and releases](#cross-compilation-and-releases)). Linux and Windows packages build in the container; electron-builder downloads Electron, the AppImage tools, and 7-Zip into `.build/desktop/cache/`. macOS packages build only on macOS, because the dmg tools and `codesign` exist only there. `dev/desktop/package --help` describes the steps: `dev/desktop/stage` puts the bundles, a minimal `package.json`, and the addon into `.build/desktop/<target>/app`, electron-builder packages that directory with `packages/app-desktop/electron-builder.config.ts`, and `dev/desktop/check` checks the Electron fuses and that the addon sits outside `app.asar`. `dev/desktop/start-test <target>` starts a package and fails when it exits early or logs an error.

The packages set the Electron [fuses](https://www.electronjs.org/docs/latest/tutorial/fuses) that turn off running the app as plain Node.js (`RunAsNode`, `NODE_OPTIONS`, `--inspect`) and load the app only from its integrity-checked `app.asar`. They are not signed with a developer certificate: on macOS they carry an ad-hoc signature, and macOS (Gatekeeper) and Windows (SmartScreen) ask the user to confirm the first start. The app has the default Electron icon.

On Ubuntu 24.04 and newer, AppArmor blocks the unprivileged user namespaces that the Chromium sandbox needs. The AppImage cannot install an AppArmor profile, so its start script detects the block and starts the app with `--no-sandbox`: on those systems the app runs without the Chromium sandbox. The AppImage uses the static AppImage runtime, which needs only `fusermount3` (installed by default on Ubuntu) and not `libfuse2`.

## Generated files

The JSON schemas, the OpenAPI document, its TypeScript client, the TypeScript declarations of the Node addon and of the WebAssembly back end, the CLI reference page of the docs site, the compiled UI messages, and the index of the example data `data/index.json` are generated and committed. `just gen` regenerates them, and `just generated-check`, part of `check-all`, fails when a committed copy is stale. Never edit them by hand.

`just data-reindex` (group `data-index` of `dev/generated`) writes `data/index.json` from the files that `git ls-files data` lists, through the `generate-data-index` binary of `packages/app-cli`; run it after adding, removing, or changing a file in `data/`, because `just generated-check` fails on a stale index and the apps offer only the files that the index lists.

`just fixtures-check` checks the reference and golden-master test fixtures against `dev/registry/reference-files.toml`.

## Testing against the reference

`dev/smoke` (host, needs Docker and `zstd`) runs the CLI over a matrix of commands, datasets and flag variants, stores the outputs in a snapshot named after the source tree that built the binary, and compares them byte for byte with a baseline snapshot, by default the one of the tip of the base branch `rust`. It reports crashes, timeouts, missing outputs and changed outputs, with the command, stderr and the changed values of each case. It also lists the cases that pass in both snapshots but print a different stderr, such as a new or repeated warning; these differences are information and never fail a case. The cases are defined in `dev/smoke.toml`, one row per command variant with its flags, datasets and expected failures; `./dev/smoke --overview` prints the cases per command and variant. `./dev/smoke --help` describes every option, the snapshot layout and the exit codes.

- `just smoke` compares the quick tier with `rust`
- `just smoke-all` compares every case with `rust`
- `just smoke-run` runs the cases without a baseline: crash, timeout and output checks only
- `just smoke-failed` runs again the cases that did not pass in the last run
- `just smoke-prune` deletes the dirty snapshots of the checkout other than the current one, and old snapshots that no branch needs. Every run also deletes the earlier current-tree snapshots of its checkout (see [Snapshot store](#snapshot-store))

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
- `--check-determinism`: run each case a second time with `--jobs=1` and mark the case nondeterministic when the outputs differ. Only the analysis commands and `pipeline` take `--jobs`; the cases of other commands run a second time unchanged
- `--jobs N`, `--timeout-scale X`: parallel cases, and a factor on every case timeout
- `--memory-budget GIB`: memory the running cases may use together, 75% of the available memory by default. Every treetime process is capped at the budget, so a case that needs more fails with an allocation error instead of exhausting the machine
- `--min-free GIB`, `--max-file-size GIB`: the disk limits (see [Disk limits](#disk-limits))

Tiers:

- `quick` (default): datasets of at most 100 sequences, the help texts and a few larger timetree cases, about 6 minutes
- `full`: every case, about 40 minutes

### Snapshot store

Without settings, each checkout keeps its own store: `snapshots/` at the checkout root (git-ignored), with the temporary baseline worktrees and the lock files under `tmp/smoke/`. With `SMOKE_STORE` set to an absolute directory in `.env` (see [Machine settings](#machine-settings)), all checkouts on the machine share one store with `snapshots/`, `worktrees/` and `locks/` in that directory. A baseline is then built once per machine, and every worktree reuses its binary and case results.

- **Snapshot id**: the first 10 hex digits of the git tree id of the commit, or `<tree>+dirty-<hash>` for a working tree with uncommitted or untracked changes. The id depends on content only, so commits with the same tree, such as a branch and the target it was fast-forwarded into, share one snapshot. `--against` resolves a git ref to the tree of its commit, and also accepts an existing snapshot id
- **Binaries**: the current binary is built in the checkout unless its snapshot already has one. A git ref baseline is built in a temporary worktree of the store with a copy of the checkout's `.env`, so it uses the same compiler cache. `dev/smoke` also sets `TREETIME_CARGO_HOME` of `dev/docker/run` to the cargo home of the checkout, so the baseline build reuses the crates that the checkout downloaded; a baseline revision whose `dev/docker/run` predates the setting downloads them again
- **Reuse**: a stored case is reused when its command, expectation, declared outputs, working directory, input files, file size limit and output storage layout are unchanged. A snapshot whose binary changed loses its stored cases
- **Outputs**: an output file of at least 1 MiB is stored compressed as `<name>.zst`; read it with `zstd -dc`. `outputs.json` holds the SHA-256 of every output after removing the volatile fields, so a comparison decompresses a file only when its hash differs
- **Locks**: a run holds an exclusive lock on a snapshot while it builds and runs it, and a shared lock while it reads the results for its report. A second run that needs the same snapshot logs the pid and checkout of the holder, waits, and then reuses the stored results. Runs on different snapshots never wait for each other. The lock is released when its holder exits, also on a crash, and an interrupted build leaves no binary behind
- **Pruning after every run**: each edit or commit gives the current tree a new snapshot that one comparison uses. After every run, `dev/smoke` deletes the snapshots that the checkout ran before as its current tree, except the current one, the baseline of the run, the snapshots of local branch tips, and snapshots used as a baseline
- **Pruning on request**: `--prune` deletes the dirty snapshots that the current checkout wrote, other than its current one, and the snapshots of committed trees older than 7 days whose tree is not the tip of any local branch. Both kinds of pruning skip snapshots that another run holds, so one worktree never deletes the work of another

Layout of one snapshot:

```text
snapshots/<id>/
  snapshot.json             commit, tree, dirty hash, ref, checkout, binary SHA-256, creation time
  bin/treetime              the binary that produced the snapshot
  last-run.json             statuses of the last run and its checkout (for --rerun-failed)
  cases/<case-id>/
    out/                    outputs; a file of at least 1 MiB is stored compressed as <name>.zst
    outputs.json            SHA-256 of every output after removing the volatile fields, and whether it is compressed
    stdout.txt, stderr.txt
    cmd.sh                  the exact command, runnable from the checkout root in the container
    result.json             status, exit code, wall and CPU time, peak memory, case key
  compare/<baseline-id>/    report.md and report.tsv of a comparison
  report.md, report.tsv     report of a --no-compare run
```

`report.md` groups failing cases by likely cause and lists value-level differences. Before the comparison, only the volatile fields `meta.updated` of `*.auspice.json` and `generated_by.version` of `*.node-data.json` are removed; numeric differences are reported, never tolerated. A `changed` case, a `regressed` case, a `still-failing` case without a declared expectation, and an `unexpected-pass` case fail the run.

Exit codes: 0 when every selected case passed, 1 when a case failed or changed, 2 when the script could not complete (build failure, git error, invalid arguments, container failure, free disk space below the limit).

### Disk limits

A full run writes large outputs; the `.dot` files of the larger datasets reach about 1 GB each. Two limits keep a run from filling the disk:

- **Free disk space** (`--min-free GIB`, default `SMOKE_MIN_FREE_GIB` in `.env`, else 20): no build starts while the file system of the checkout or the store has less free space. While cases run, the free space of the store is checked every second; below the limit, the running cases are stopped, their partial outputs deleted, no further case starts, and the run exits with code 2. The cases that completed keep their results, so the next run continues from them
- **File size** (`--max-file-size GIB`, default 8): every treetime process may write files up to this size. A larger file stops the process with `SIGXFSZ`, and the case fails with that reason

### TreeTime v0 and Python scripts

TreeTime v0 (`packages/legacy/treetime`), the golden-master capture scripts, and the research scripts in `test_scripts/` share one Python environment, defined in `pixi.toml` and locked in `pixi.lock`. It pins every package and uses the OpenBLAS 0.3.34 release of conda-forge on a single thread. The Rust build links OpenBLAS develop `de5d4914`, 432 commits after 0.3.34, with LAPACK compiled by flang, so v0 and v1 run close but not identical linear algebra kernels. The environment puts v0 of the checkout on `PYTHONPATH`, so both the v0 `treetime` command and `import treetime` in scripts work.

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

- `dev` (`just build`, `just run`, the tests): unoptimized workspace crates, dependencies at `opt-level = 2`, full debug info. Every build for `x86_64-unknown-linux-gnu`, the development builds included, compiles Rust and C code for the CPU level of the shipped x86_64 binaries, x86-64-v2 (`.cargo/config.toml`)
- `dev-opt` (`just run-dev-opt`): `dev` with optimized workspace crates, for long runs on real datasets; rebuilds are slower
- `release` (`just build-release`, `just run-release`, `just example`, `just smoke`): optimized and fast to rebuild, without LTO
- `dist` (`just build-dist`, `just cross`, `just build-napi dist`, `just cross-napi`, nightly releases): the shipped binary and Node addon, with fat LTO, one codegen unit, and the CPU and link flags of the shipped builds from `export_dist_flags` of `dev/lib/dist-flags.sh` (`-C target-cpu=x86-64-v2` on x86-64, no frame pointers). These flags replace the target rustflags of the Cargo config files, so the native dist builds link with rust-lld. The glibc binaries of `build-dist`, `build-profiling`, `build-profiling-heap`, and `bench` link the GCC unwinder statically and need no `libgcc_s.so.1`, as the cross builds do; the Node addon keeps the shared unwinder of the Electron process
- `profiling` (`just build-profiling`, `just profile`): `dist` with full debug info; without frame pointers, profilers unwind call stacks with DWARF (`dev/profile` does) or LBR
- `profiling-heap` (`just build-profiling-heap`): `profiling` with jemalloc heap profiling compiled in, for heap profiles only
- `bench` (`just bench`): the `dist` settings

## Performance

- `just bench` runs the benchmarks
- `just build-dist treetime` builds the shipped binary into `.out/treetime`; use it with `hyperfine`, which the container provides
- `just profile treetime -- <args>` (host) samples a profile with samply or perf; read `dev/profile --help` first

### Allocator

On x86_64 and aarch64 Linux, the CLI and the API server allocate through jemalloc. jemalloc is built without a symbol prefix (the `override_allocator_on_supported_platforms` feature of `tikv-jemalloc-sys`), so the binary defines `malloc`, `free`, and the other allocator functions of the C library itself. Rust code, the C libraries (OpenBLAS, liblzma, libzstd), and the C library all allocate through jemalloc. The glibc binaries export these functions, so glibc binds its own calls to them. The static musl binaries keep the musl allocator only for the internal allocations that musl also frees itself. Other targets and the Node addon use the allocator of the system.

- `MALLOC_CONF` configures jemalloc, for example `MALLOC_CONF=dirty_decay_ms:1000`. The settings of the glibc allocator (`MALLOC_ARENA_MAX`, the `glibc.malloc` tunables) have no effect
- `dev/cross/check` fails when a Linux binary does not take these functions from jemalloc

### Heap profiles

Heap profilers that replace `malloc` through `LD_PRELOAD` (heaptrack) do not see the allocations of the CLI, because the binary binds its calls to its own `malloc`. jemalloc's own heap profiler sees them. It is compiled into the `profiling-heap` build only, because it changes the allocator code: take CPU profiles and timings from the `profiling` and `dist` builds.

```bash
./dev/docker/run just build-profiling-heap
mkdir -p tmp/heap
MALLOC_CONF=prof:true,prof_gdump:true,prof_final:true,lg_prof_sample:17,prof_prefix:tmp/heap/treetime \
  .build/container/cargo/profiling-heap/treetime <command> <args>
peak="$(ls -v tmp/heap/*.u*.heap | tail -n 1)"
jeprof --collapsed .build/container/cargo/profiling-heap/treetime "${peak}" > tmp/heap/peak.folded
flamegraph.pl --countname=bytes tmp/heap/peak.folded > tmp/heap/peak.svg
```

- `prof_gdump` writes a `.u<n>.heap` dump each time the heap grows past its previous maximum, so the dump with the highest number is the closest to the peak. A run on a few hundred sequences writes thousands of dumps and several hundred MB; delete the others once the peak dump is chosen. `prof_final` writes one more dump (`.f.heap`) at exit
- `lg_prof_sample:17` samples once per 128 KiB allocated on average (the default is 512 KiB)
- `--collapsed` keeps full Rust paths and inlined frames. `--text` shortens Rust names to their last segment and charges every byte to jemalloc's sampler, which makes it hard to read
- `jeprof --base=<earlier>.heap` subtracts an earlier dump, and `--alloc_space` reports total allocation volume when the run had `prof_accum:true`
- `jeprof` comes with the `libjemalloc-dev` package on Debian and Ubuntu; `flamegraph.pl` is in [FlameGraph](https://github.com/brendangregg/FlameGraph)

## Dependencies

Every dependency release must be at least seven days old before the project adopts it:

- Bun refuses younger npm packages (`minimumReleaseAge` in `bunfig.toml`)
- Cargo has the age check as an unstable feature until Rust 1.100. `just deps-update` and `just deps-upgrade` enable it for their resolution, and `just deps-age` fails on a lockfile entry younger than seven days
- Python packages: `exclude-newer` in `pixi.toml` makes `pixi lock` skip younger releases
- Tools in `.config/mise.toml`: `minimum_release_age` makes `just tools-outdated` list only newer releases that are at least seven days old. mise does not filter an exact version, so check the date of a version you type by hand
- Base images in `dev/docker/` follow the same rule by hand
- Releases of the [binarylandia](https://github.com/binarylandia) GitHub organization (the prebuilt OpenBLAS, liblzma, and libzstd archives, the compilers, the patched `wasm-bindgen` CLI) are the project's own builds and need no waiting time; `minimum_release_age_excludes` in `.config/mise.toml` lists their repositories. The sources that such a release builds must still be at least seven days old
- mise itself: the images and CI install the version in `dev/docker/files/mise-version`, verified by its line in `dev/docker/files/checksums`. `min_version` in `.config/mise.toml` is the oldest mise that reads the configuration; raise it when the configuration needs a newer mise
- After changing a tool in `.config/mise.toml`, `just tools-lock <tool>` locks only that tool, and after changing an entry of the images in `.config/mise.image.toml`, `just tools-lock-image <tool>` locks it in `.config/mise.image.lock`. Locking calls the GitHub API, which allows 60 anonymous requests per hour. With `MISE_GITHUB_TOKEN` set, mise authenticates instead, on the host and in the container alike (`dev/docker/run` forwards it), for example `MISE_GITHUB_TOKEN="$(gh auth token)" ./dev/docker/run just tools-lock <tool>`

Rust dependencies are pinned exactly in the workspace `Cargo.toml`, JavaScript dependencies exactly in the manifests, with shared packages in the Bun catalog of the root `package.json`. The React packages stay on 19.x, because Auspice runs in-process and each React major needs its own check of the Auspice modules TreeTime imports; `just lint-ts` enforces this.

`just deps-update` changes `Cargo.lock` in any checkout: after a new dependency in `Cargo.toml`, `just deps-update --workspace` adds it to the lockfile without updating the other entries. The upgrade recipes `just deps-upgrade` and `just deps-upgrade-ts` run in the main checkout only. `just audit` checks both dependency graphs against the security advisory databases.

## Cross-compilation and releases

`just cross` (host) builds the shipped CLI (`dist` profile) for every target in `dev/cross/targets` in its cross image, in parallel, into `.out/`. The targets divide the CPU budget `TREETIME_DOCKER_CPUS` between them (see [Machine settings](#machine-settings)). `just cross --target=aarch64-apple-darwin` builds one target. Without `just` on the host, run `./dev/cross/all treetime`.

`just cross-napi` (host) builds the Node addon of the desktop app the same way, with the same profile and CPU flags, for every target in `dev/cross/targets-desktop`, the targets Electron supports (no musl), into `.out/app-napi-<target>.node` (`./dev/cross/all --lib app-napi` without `just`). The Windows addon is built with llvm-mingw, like the CLI, and loads into Electron, which is built with MSVC: Node-API reaches the addon through the functions that `napi-sys` looks up in the running executable, not through an import library.

Each target family has its own cross image. mise installs every toolchain component and library archive of the images from `.config/mise.image.toml` at the URL and checksum of `.config/mise.image.lock` (`dev/docker/files/install-image-tools`), and the libraries link through pkg-config, as in the native image (see [Build container](#build-container-recommended)):

- Linux (`dev/docker/cross-linux.dockerfile`): clang and lld of LLVM 22 ([binarylandia/build_llvm](https://github.com/binarylandia/build_llvm)) on the sysroot, crt files, libgcc, and libstdc++ of a GCC 16 toolchain of [binarylandia/build_crosstool-ng](https://github.com/binarylandia/build_crosstool-ng): glibc 2.17 for the gnu targets, musl 1.2.6 for the musl targets. The musl binaries are static
- Windows (`dev/docker/cross-windows.dockerfile`): llvm-mingw of the same LLVM release, with the static mingw-w64 sysroot on the Universal C Runtime (UCRT), compiler-rt, and libunwind, for `x86_64-pc-windows-gnullvm` and `aarch64-pc-windows-gnullvm`, and the OpenBLAS, liblzma, and libzstd archives built with it. The binaries run on Windows 10 and later
- Windows 7 (`dev/docker/cross-windows7.dockerfile`, target `x86_64-win7-windows-gnu`, not in `dev/cross/targets`): the GCC 16 MinGW toolchain of build_crosstool-ng on the msvcrt C runtime, which Windows 7 has without updates, the msvcrt OpenBLAS of build_openblas, and the liblzma and libzstd archives of staticbourg built with this GCC toolchain. Rust ships no standard library for this tier 3 target, so `dev/cross/build` compiles it from the `rust-src` component with `-Zbuild-std` (and `RUSTC_BOOTSTRAP=1`, which the stable toolchain needs for that flag). The CLI only: Electron does not run on Windows 7. Build it with `just cross --target=x86_64-win7-windows-gnu`
- WebAssembly (`dev/docker/cross-wasm.dockerfile`, target `wasm32-unknown-emscripten`, not in `dev/cross/targets`): Emscripten, the Rust toolchain with the nightly of the threaded module, the WebAssembly archives of OpenBLAS, liblzma, and libzstd, and the tools of `just build-web2`, without the native toolchains, the lint toolchains, and the other development tools of the native image. The GCC of Debian compiles the build scripts. The web deploy builds in it, and so can you: `CROSS_COMPILE=wasm32-unknown-emscripten ./dev/docker/run just build-web2`. The development image keeps Emscripten, for the WebAssembly lints and tests of `just check`
- macOS (`dev/docker/cross-darwin.dockerfile`): the llvm flavor of osxcross ([binarylandia/build_osxcross](https://github.com/binarylandia/build_osxcross)) with the macOS 11.1 SDK, which runs clang and `ld64.lld` of the LLVM release, and compiler-rt. The deployment targets are macOS 10.12 (x86_64) and 11.0 (arm64)

`dev/cross/check` fails when a binary or addon needs a library that is not part of the operating system: anything beyond glibc 2.17 on the Linux gnu targets (the addon may also need `libgcc_s.so.1`, as Electron does), any library on the musl targets, a DLL other than the Windows system DLLs (`KERNEL32.dll`, `api-ms-win-crt-*.dll`, and others) on Windows, a DLL other than the Windows 7 system DLLs and `msvcrt.dll`, or a function that Windows 8 or later added, on Windows 7, and a library outside `/usr/lib` and `/System/Library` on macOS. It also fails when a macOS binary names another minimum macOS version than its deployment target. `dev/cross/all --run` runs each Linux binary under QEMU and each x86_64 Windows binary under Wine; macOS and aarch64 Windows binaries are only built and checked there. For each Linux binary it also runs `dev/cross/check-cpus`, which runs the binary under QEMU on the oldest CPU of its baseline (Nehalem, Cortex-A53) with every compression format, and on x86_64 below the baseline, where the binary must stop with its CPU message. CI runs `dev/cross/check` after each build and `dev/cross/check-cpus` in the `cpu-baseline` jobs, and it runs the Linux binaries on distribution images for x86_64 and arm64 (`compat-linux`), the Windows binaries on `windows-2022` and `windows-11-arm` (`compat-windows`), and the macOS binaries on Intel and Apple Silicon macOS runners (`compat-macos`).

The release binaries run on these CPUs; on an older x86_64 CPU they stop with a message that names the missing level. The prebuilt archives they link follow the same levels: build_openblas and staticbourg compile them with the C CPU flags of `export_dist_flags`, which the image build checks in the build record of every archive, and OpenBLAS also selects its kernels at run time.

- x86_64 (every operating system, Windows 7 included): x86-64-v2 (SSE4.2, POPCNT, CMPXCHG16B; Intel Nehalem, AMD Bulldozer, or newer)
- Linux aarch64: Armv8.0, such as the Cortex-A53 and Cortex-A72 of the Raspberry Pi 3 and 4. Atomic operations go through outline atomics, which use LSE instructions when the CPU has them
- Windows aarch64: the LSE atomics of Armv8.1, which every device that runs Windows 11 on Arm has
- macOS aarch64: Apple M1 or newer

### Nightly releases

Prerelease builds are published in [neherlab/treetime-nightly](https://github.com/neherlab/treetime-nightly/releases). A nightly publishes every CLI binary and desktop package that builds, without the linkage checks and start tests of the CI workflow.

1. `.github/workflows/nightly.yml` on `rust` runs on every push to `rust` and releases its tip. Runs do not overlap: a push during a running nightly waits, and only the newest waiting run starts, so one release can cover several commits. `.github/workflows/schedule-nightly.yml` on `master` also calls `nightly.yml` daily at 04:00 UTC, which skips when the latest nightly already has the tip of `rust`, so it releases only a commit whose push run failed
2. `nightly.yml` builds the cross-compilation matrix of `.github/workflows/cli-build.yml`, and the desktop apps with `.github/workflows/desktop-build.yml`: the addon of each target in its cross image, the bundles, and the packages. A failed target does not block the others
3. `nightly.yml` writes `SHA256SUMS` over every asset of the release and attests these checksums with `actions/attest`, which signs with the identity of the workflow run. The attestation needs the `id-token: write` and `attestations: write` permissions. A called workflow can only lower the permissions of its caller, so the `nightly` job of `schedule-nightly.yml` grants them too
4. `dev/publish-nightly` creates the prerelease with the built binaries, the desktop packages, `examples.zip`, which holds the files that `data/index.json` lists and the index itself, and `SHA256SUMS`. Its notes list the changelog since the previous nightly, followed by a downloads table that links every asset, one row per platform, and the command that verifies a download: `gh attestation verify <file> --repo neherlab/treetime`

`./dev/trigger-nightly` starts a nightly by hand. The version format is `<cargo-version>-nightly.<YYYYMMDD>T<HHMMSS>Z+<short-sha>`.

`./dev/trigger-nightly --only cli` runs only the CLI release, `--only desktop` only the desktop release; push runs and the scheduled run do both. A desktop-only release has no CLI binaries, so the next scheduled run still builds the CLI for the same commit, and a desktop-only run also runs when `rust` is unchanged.

### Web app deploys

`.github/workflows/web.yml` builds the web2 app (`just build-web2`) on every push in the WebAssembly image (`CROSS_COMPILE=wasm32-unknown-emscripten`, see [Cross-compilation and releases](#cross-compilation-and-releases)) and uploads it with `dev/publish-web` to the preview host, a CloudFront distribution in front of an S3 bucket in the AWS account `preview-deployments` of neherlab. `just build-web2` compiles the single-threaded and the threaded module at the same time. `dev/publish-web` uploads the build to the first site and copies it to the other sites inside the bucket. Each upload goes to several sites:

- `https://neherlab--treetime--commit-<short-sha>.previews.neherlab.click`: the commit
- `https://neherlab--treetime--branch-<branch>.previews.neherlab.click`: the latest upload of the branch
- `https://neherlab--treetime--pr-<number>.previews.neherlab.click`: the latest upload of each open pull request of the branch; the run adds or updates a comment with the link in the pull request
- `https://treetime.neherlab.click`: the latest push to `rust`

Opening or reopening a pull request also starts a run, so a pull request gets its site without a new push. Pull requests from forks build without an upload, because workflows of forks get no secrets. The upload reads the access key, the region, the bucket, and the distribution from the organization secrets `PREVIEWS_NEHERLAB_CLICK_AWS_ACCESS_KEY_ID`, `PREVIEWS_NEHERLAB_CLICK_AWS_SECRET_ACCESS_KEY`, `PREVIEWS_NEHERLAB_CLICK_AWS_DEFAULT_REGION`, `PREVIEWS_NEHERLAB_CLICK_AWS_S3_BUCKET`, and `PREVIEWS_NEHERLAB_CLICK_AWS_CLOUDFRONT_DISTRIBUTION_ID`. `dev/publish-web --help` describes the stored files, their caching, and what the edge functions of the distribution do; the edge functions themselves are not in this repository.

### Docs deploys

`.github/workflows/docs.yml` builds the docs site (`just build-docs`) on every push and uploads it with the same composite action as the web app (`.github/actions/publish-preview`). Each docs build is a version folder of its site, today only `dev`, the tip of `rust`:

- `https://neherlab--treetime-docs--commit-<short-sha>.previews.neherlab.click/dev`, `...--branch-<branch>...`, and `...--pr-<number>...` (with a pull request comment marked `<!-- docs-preview -->`)
- `https://docs.treetime.neherlab.click/dev`: the latest push to `rust`, from the bucket folder `treetime-docs`

Each upload also stores the site root files of `dist-root/`: `versions.json`, the list of versions, and an `index` page that leads to the preferred version. The version select and the banner of every build read `/versions.json` at run time, so an older build learns about a newer version without a rebuild. `packages/app-docs/versions.json` is the source of the list; `TREETIME_DOCS_VERSION` names the version of a build (default `dev`), and the build fails when the list has no entry for it. Adding a release version means building the docs at the release tag with `TREETIME_DOCS_VERSION=<version>`, uploading with `dev/publish-web --prefix <version>`, and adding the entry, marked `preferred`, to `versions.json`.

## Continuous integration

`.github/workflows/cli.yml` runs on pull requests and on pushes to `rust`:

- The check groups of `just check-all` in parallel jobs, except `hawk`, which needs a second full compilation; run it locally
- On pushes to `rust`: the release builds of every target and compatibility runs on Linux distributions, macOS, and Windows

`.github/workflows/web.yml` builds and uploads the web2 app on every push (see Web app deploys), and `.github/workflows/docs.yml` the docs site (see Docs deploys).

The CI jobs pull the container images from Docker Hub by the hash of their build inputs, and build them when the inputs changed. Only pushes to `rust` publish images, with zstd-compressed layers, which compress and unpack faster than gzip and need Docker 23 or later to pull. The layers of each image go from the inputs that change least to those that change most, so a tool update in `.config/mise.toml` rebuilds only the last layers.
