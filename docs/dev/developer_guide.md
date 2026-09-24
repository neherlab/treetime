# Developer guide

This guide describes how to set up a development environment, build and test TreeTime, and maintain the project. It assumes basic familiarity with TreeTime from a user perspective.

The `justfile` is the entry point for every routine task. `just` lists the recipes by group, and each recipe has a one-line description.

## Setup

TreeTime builds in two ways: in the build container (recommended), or directly on the host. Both run the same `just` recipes with the same tool versions, which `mise.toml` pins and `mise.lock` locks by URL and checksum.

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

- `KACHE_STORE`: directory of the [kache](https://github.com/kunobi-ninja/kache) compiler cache store. Builds and clippy then compile through kache, which shares compiled crates across worktrees
- `TREETIME_PORTLESS`: use of the [portless](https://github.com/vercel-labs/portless) proxy by the web dev server (see [Web app](#web-app))

## Everyday commands

| Task                                             | Command                               |
| ------------------------------------------------ | ------------------------------------- |
| List the recipes                                 | `just`                                |
| Fast checks (format, clippy, TypeScript)         | `just check`                          |
| Every check; must pass before merging            | `just check-all`                      |
| Fast lint fixes and format (stage changes first) | `just fix`                            |
| Every lint fix, dylint included, and format      | `just fix-all`                        |
| Build                                            | `just build` (`just b`)               |
| Run the CLI                                      | `just run treetime ancestral --help`  |
| Rust tests, optionally filtered                  | `just test-rs [filter]` (`just t`)    |
| TypeScript tests                                 | `just test-ts`                        |
| Clippy                                           | `just lint-rs` (`just l`)             |
| TypeScript lints                                 | `just lint-ts`                        |
| Format one toolchain                             | `just fmt-rs`, `fmt-ts`, `fmt-other`  |

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
kb/                  Knowledge base
data/                Example datasets
```

## Apps

### Web app

`just up` starts the API server and the Vite dev server in the foreground until Ctrl-C. `just health` and `just status` report whether they run and which commit they were started from. Each checkout resolves its own ports, so worktrees run side by side; `TREETIME_API_PORT` and `TREETIME_WEB_PORT` override them.

When the portless proxy runs on the machine, the web server also registers `https://treetime.localhost` in the main checkout and `https://<branch>.treetime.localhost` in a linked worktree.

### Desktop app

`just desktop` starts the Electron app in development mode. In the container it needs the host display: `TREETIME_DOCKER_X11=1 ./dev/docker/run just desktop`.

## Generated files

The JSON schemas, the OpenAPI document, its TypeScript client, and the CLI reference documentation are generated and committed. `just gen` regenerates them, and `just generated-check`, part of `check-all`, fails when a committed copy is stale. Never edit them by hand.

`just fixtures-check` checks the reference and golden-master test fixtures against `dev/registry/reference-files.toml`.

## Testing against the reference

`dev/smoke` (host, needs Docker) runs the CLI over a matrix of commands, datasets and flag variants, stores the outputs in a snapshot under `snapshots/<id>/` named after the git state that built the binary, and compares them byte for byte with a baseline snapshot, by default the one of the `rust` branch. It reports crashes, timeouts, missing outputs and changed outputs, and writes `report.md` and `report.tsv` into `snapshots/<id>/compare/<baseline-id>/`. The cases are defined in `dev/smoke.toml`, one row per command variant with its flags, datasets and expected failures; `./dev/smoke --overview` prints the cases per command and variant. `./dev/smoke --help` describes the options, the snapshot layout and the exit codes.

- `just smoke` compares the quick tier (datasets of at most 100 sequences) with `rust`
- `just smoke-all` compares every case with `rust`
- `just smoke-run` runs the cases without a baseline: crash, timeout and output checks only
- `just smoke-failed` runs again the cases that did not pass in the last run
- `just smoke-prune` deletes the snapshots of dirty working trees other than the current one

```bash
./dev/smoke                    # quick tier against rust
./dev/smoke --tier full        # every case
./dev/smoke --against <ref>    # another commit, branch or snapshot id as the baseline
./dev/smoke --no-compare       # statuses only
./dev/smoke --only <regex>     # cases whose id matches
./dev/smoke --rerun-failed     # cases that did not pass in the last run
```

- `dev/docker/python treetime ...` runs TreeTime v0 from `packages/legacy`

## Performance

- `just bench` runs the benchmarks
- `just build-release treetime` builds the release binary into `.out/treetime`; use it with `hyperfine`, which the container provides
- `just profile treetime -- <args>` (host) samples a profile with samply or perf; read `dev/profile --help` first

## Dependencies

Every dependency release must be at least seven days old before the project adopts it:

- Bun refuses younger npm packages (`minimumReleaseAge` in `bunfig.toml`)
- Cargo has the age check as an unstable feature until Rust 1.100. `just deps-update` and `just deps-upgrade` enable it for their resolution, and `just deps-age` fails on a lockfile entry younger than seven days
- Tools in `mise.toml` and base images in `dev/docker/` follow the same rule by hand; `just tools-outdated` lists newer tool releases

Rust dependencies are pinned exactly in the workspace `Cargo.toml`, JavaScript dependencies exactly in the manifests, with shared packages in the Bun catalog of the root `package.json`. The React packages stay on 18.x, because Auspice runs in-process and requires it; `just lint-ts` enforces this.

The dependency recipes run in the main checkout only. `just audit` checks both dependency graphs against the security advisory databases.

## Cross-compilation and releases

`just cross` (host) builds the release CLI for every target in `dev/cross/targets` in its cross image, in parallel, into `.out/`. `just cross --target=aarch64-apple-darwin` builds one target. Without `just` on the host, run `./dev/cross/all treetime`.

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
