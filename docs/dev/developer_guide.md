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

Host builds go to `.build/host/`. The dylint and IQ-TREE tools are available on Linux only; run the recipes that need them in the container on other hosts.

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

The interface is built from [shadcn/ui](https://ui.shadcn.com) components in the `base-vega` style, which run on Base UI. The components are vendored into `packages/app-ui/src/ui/`, one file per component, and the app composes them; it does not style raw elements. The theme tokens (colors, radius, fonts) live in `packages/app-ui/src/theme.css`, and `packages/app-ui/src/ui/shadcn.css` is the stylesheet of the `shadcn` npm package that defines the variants the components use.

To add or update a component, fetch its source from the registry, `https://ui.shadcn.com/r/styles/base-vega/<name>.json` (the `files[0].content` field), and apply the rewrites the shadcn CLI would apply:

- **Imports**: `@/registry/base-vega/ui/<name>` becomes `./<name>`, and `cn` comes from `./cn`
- **Icons**: each `IconPlaceholder` element becomes the `lucide-react` icon named in its `lucide` attribute
- **Classes**: `cn-font-heading` becomes `font-heading`, and the other `cn-*` marker classes are removed
- **Comments**: removed, as in all TypeScript source

The vendored files keep the upstream code shape, so `oxlint.config.ts` turns off the style rules for `packages/app-ui/src/ui/*.tsx` while the Tailwind class check stays on.

Library hooks cover the behavior around the components: `@tanstack/react-hotkeys` for keyboard shortcuts, `react-dropzone` for file drops, `use-stick-to-bottom` for the following log, `@tanstack/react-table` for sortable tables, `@tanstack/react-pacer` for debouncing, `cmdk` for the command palette, and `@mantine/hooks` for the clipboard, element sizes, media queries, and the file dialog.

### Web app

`just up` starts the API server and the Vite dev server in the foreground until Ctrl-C. `just health` and `just status` report whether they run and which commit they were started from.

`just serve` builds the web app and the API server for production and runs the server in the foreground until Ctrl-C. The server serves the built web app (`packages/app-web/dist`, through the `STATIC_DIR` environment variable) and the API from one port.

Each checkout resolves its own ports, so worktrees run side by side; `TREETIME_API_PORT`, `TREETIME_WEB_PORT`, and `TREETIME_SERVE_PORT` override them.

#### Storage

The web app has no database. The server keeps each run in its own folder `<runs-dir>/<run-id>/`: the run record `run.json`, the event log `events.jsonl`, the uploaded files in `inputs/`, and the results in `out/`. The browser keeps only the unfinished form, in `localStorage`. The server reads example datasets from `--data-dir`.

Each app mode keeps its runs in its own directory:

- `tmp/app/web-dev/runs`: development servers (`just up`)
- `tmp/app/web-prod/runs`: production server (`just serve`)
- `tmp/app/desktop-dev/runs`: desktop app in development mode (`just desktop`), with its diagnostics in `tmp/app/desktop-dev/diagnostics`
- `tmp/app/desktop-prod/runs`: production build of the desktop app started from the checkout (`just desktop-prod`), with its diagnostics in `tmp/app/desktop-prod/diagnostics`
- `<userData>/runs`: installed desktop app, in the Electron user data directory of the platform

A deployment passes `--data-dir` and `--runs-dir` to `treetime-server` itself.

When the portless proxy runs on the machine, the web server also registers `https://treetime.localhost` in the main checkout and `https://<branch>.treetime.localhost` in a linked worktree.

#### Event streams and HTTP/2

The web app follows changes through server-sent event streams: `GET /api/events` for changes to any run, and `GET /api/runs/{id}/events` for each run it follows. Every open stream holds one connection. Over HTTP/1.1, browsers allow at most 6 connections per host, shared by all tabs, so more open streams stall every further request to the server. Over HTTP/2, all requests share one connection.

The API server speaks HTTP/1.1 without TLS, and browsers use HTTP/2 only over TLS. A deployment therefore puts a reverse proxy in front of the server that terminates TLS and serves HTTP/2 to the browser, for example nginx with `listen 443 ssl; http2 on;`. In development, the portless proxy (`https://treetime.localhost`) does this by default. Without such a proxy, at most 6 streams stay open per host across all tabs.

### Desktop app

`just desktop` starts the Electron app in development mode. `just desktop-prod` builds the desktop app for production (the release Node addon and the bundled renderer) and starts that build from the checkout. In the container both need the host display: `TREETIME_DOCKER_X11=1 ./dev/docker/run just desktop`.

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
- `dist` (`just build-dist`, `just cross`, nightly releases): the shipped binary, with fat LTO and one codegen unit
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
- Tools in `.config/mise.toml` and base images in `dev/docker/` follow the same rule by hand; `just tools-outdated` lists newer tool releases
- After changing a tool in `.config/mise.toml`, `just tools-lock <tool>` locks only that tool. Locking calls the GitHub API, which allows 60 anonymous requests per hour. With `MISE_GITHUB_TOKEN` set, mise authenticates instead, on the host and in the container alike (`dev/docker/run` forwards it), for example `MISE_GITHUB_TOKEN="$(gh auth token)" ./dev/docker/run just tools-lock <tool>`

Rust dependencies are pinned exactly in the workspace `Cargo.toml`, JavaScript dependencies exactly in the manifests, with shared packages in the Bun catalog of the root `package.json`. The React packages stay on 19.x, because Auspice runs in-process and each React major needs its own check of the Auspice modules TreeTime imports; `just lint-ts` enforces this.

The dependency recipes run in the main checkout only. `just audit` checks both dependency graphs against the security advisory databases.

## Cross-compilation and releases

`just cross` (host) builds the shipped CLI (`dist` profile) for every target in `dev/cross/targets` in its cross image, in parallel, into `.out/`. `just cross --target=aarch64-apple-darwin` builds one target. Without `just` on the host, run `./dev/cross/all treetime`.

The release binaries require an x86_64 CPU with AVX2 (Haswell or newer) and, on Linux aarch64, ARMv8.2.

### Nightly releases

Prerelease builds are published in [neherlab/treetime-nightly](https://github.com/neherlab/treetime-nightly/releases). A nightly publishes every target that builds; the `x86_64-unknown-linux-gnu` binary is required.

1. `.github/workflows/schedule-nightly.yml` on `master` runs daily at 04:00 UTC and calls `.github/workflows/nightly.yml` on `rust`, which skips when `rust` has no new commits
2. `nightly.yml` builds the cross-compilation matrix of `.github/workflows/cli-build.yml`
3. `dev/publish-nightly` creates the prerelease with the built binaries

`./dev/trigger-nightly` starts a nightly by hand. The version format is `<cargo-version>-nightly.<YYYYMMDD>T<HHMMSS>Z+<short-sha>`.

## Continuous integration

`.github/workflows/cli.yml` runs on pull requests and on pushes to `rust`:

- The check groups of `just check-all` in parallel jobs, except `hawk`, which needs a second full compilation; run it locally
- On pushes to `rust`: the release builds of every target and compatibility runs on Linux distributions, macOS, and Windows

The CI jobs pull the container images from Docker Hub by the hash of their build inputs, and build them when the inputs changed. Only pushes to `rust` publish images.
