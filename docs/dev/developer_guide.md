# Contributing to Treetime

This guide describes how to set up a development environment, build Treetime, contribute to the codebase, and maintain the project.

> ⚠️ This is a document for developers, maintainers of the project as well as for the most curious users. It assumes basic familiarity with Treetime program from a user perspective. If you are not familiar with Treetime, read the User Documentation first.

## Setup Development Environment

### Docker (recommended)

The recommended way to develop Treetime is using the Docker-based dev environment. It provides a consistent setup with all Rust, Node.js, and build dependencies pre-installed.

Requirements: bash, Docker

```bash
# First run builds the Docker image (takes a few minutes)
./dev/docker/run ./dev/dev l

# Rebuild the image after Dockerfile changes
DOCKER_FORCE_REBUILD=1 ./dev/docker/run echo ok
```

All `./dev/dev` commands should be run inside the Docker container via `./dev/docker/run`.

### Native (alternative)

This guide assumes you are using Ubuntu 24.04, but the instructions will likely be similar for other Linux and Unix-like operating systems.

Treetime is written in Rust, and the standard `rustup` & `cargo` workflow is used.

```bash
sudo apt-get update
sudo apt-get install \
  bash \
  build-essential \
  clang \
  curl \
  gcc \
  gfortran \
  git \
  libbz2-dev \
  libclang-dev \
  liblzma-dev \
  libopenblas-dev \
  libssl-dev \
  libzstd-dev \
  make \
  pkg-config \
  protobuf-compiler \
  zlib1g-dev \

# (optional) To enable PNG image output ("png" cargo feature, see below), install:
sudo apt-get install libfontconfig1-dev

# Install Rustup, the Rust version manager (https://www.rust-lang.org/tools/install)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y

# Add Rust tools to the $PATH
export PATH="$PATH:$HOME/.cargo/bin"

# Additional tools for testing and maintenance
cargo -q install --locked cargo-nextest cargo-edit cargo-audit
```

## Obtain Source Code

Treetime is an open-source project licensed under the MIT license, and its source code is available on GitHub.

```bash
git clone https://github.com/neherlab/treetime
```

If you are a team member, use the SSH URL:

```bash
git clone git@github.com:neherlab/treetime.git
```

If you are not a team member but want to contribute, make a [fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo) and clone your forked repository instead. You can then submit changes as a [Pull Request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request).

## Project Structure

```
packages/
  app-api/           Rust service layer (progress, commands)
  app-cli/           Rust CLI binary
  app-contracts/     TypeScript bridge types
  app-desktop/       Electron shell
  app-napi/          Rust napi addon (native Node.js module)
  app-server/        Rust HTTP API server
  app-ui/            Shared React components
  app-web/           Vite SPA (browser frontend)
  legacy/            Python v0 reference implementation
  treetime/          Core Rust library
  treetime-*/        Core supporting crates
  util-phyloxml/     PhyloXML format library
  util-usher-mat/    UShER MAT format library
```

## Dev Shortcuts

All commands are run via `./dev/docker/run ./dev/dev <shortcut>`. Multiple shortcuts can be chained: `./dev/dev l f t` runs lint, format, then test.

### Rust

| Shortcut   | Command          | Description                           |
| ---------- | ---------------- | ------------------------------------- |
| `b`        | build            | Debug build                           |
| `br`       | build-release    | Release build                         |
| `r <bin>`  | run              | Run binary in debug mode              |
| `rr <bin>` | run-release      | Run binary in release mode            |
| `l`        | lint             | Clippy (includes build)               |
| `lf`       | lint-fix         | Clippy with auto-fix                  |
| `lc`       | lint-ci          | Clippy with -Dwarnings                |
| `lD`       | lint-deps        | Unused dependency check (cargo-shear) |
| `f`        | format           | rustfmt                               |
| `fc`       | format-check     | rustfmt check only                    |
| `t`        | test-all         | All tests (nextest)                   |
| `tu`       | test-unit        | Unit tests only                       |
| `ti`       | test-integration | Integration tests only                |
| `tl`       | test-list        | List all tests                        |
| `q`        | quality          | lint-fix + format + test              |
| `B`        | bench-all        | Run benchmarks                        |
| `U`        | upgrade-deps     | Upgrade Cargo dependencies            |
| `cov`      | coverage         | HTML coverage report                  |
| `w`        | why              | Dependency tree                       |
| `L`        | list             | List all binaries, examples, tests    |

### Desktop and Web Apps

| Shortcut | Command           | Description                                                     |
| -------- | ----------------- | --------------------------------------------------------------- |
| `d`      | desktop-start     | Start Electron desktop app (Vite renderer + napi addon via IPC) |
| `db`     | desktop-build     | Production build of desktop app                                 |
| `a`      | app-start         | Start web app (Bacon server auto-restart + Vite HMR)            |
| `ab`     | app-build         | Production build of web app                                     |
| `ar`     | app-release       | Release build of web app                                        |
| `arw`    | app-release-watch | Release build with watch mode                                   |

### JavaScript/TypeScript (bulk, all packages)

| Shortcut | Command         | Description                                |
| -------- | --------------- | ------------------------------------------ |
| `ji`     | js-install      | Install all npm dependencies (bun install) |
| `ju`     | js-upgrade      | Upgrade all npm dependencies to latest     |
| `jl`     | js-lint         | Lint all TypeScript (oxlint)               |
| `jlf`    | js-lint-fix     | Lint with auto-fix                         |
| `jf`     | js-format       | Format all TypeScript (prettier)           |
| `jfc`    | js-format-check | Format check only                          |
| `jc`     | js-check        | Typecheck all TypeScript (tsc --noEmit)    |
| `jt`     | js-test         | Run all TypeScript tests (vitest)          |

## Build and Run

### CLI

```bash
v="flu/h3n2/20"

# Run in debug mode
./dev/docker/run ./dev/dev r treetime -- ancestral --method-anc=parsimony --outdir="tmp/ancestral/$v" --tree="data/$v/tree.nwk" "data/$v/aln.fasta.xz"

# Run in release mode
./dev/docker/run ./dev/dev rr treetime -- ancestral --method-anc=marginal --model=jc69 --outdir="tmp/ancestral/$v" --tree="data/$v/tree.nwk" "data/$v/aln.fasta.xz"

# Clock estimation
./dev/docker/run ./dev/dev r treetime -- clock --tree="data/$v/tree.nwk" --dates="data/$v/metadata.tsv" --outdir="tmp/clock/$v"
```

> 💡 Set variable `v` to a different path in the `data/` directory to try other example inputs. List all examples with:
>
> ```
> find data -type f -name "tree.nwk" -exec dirname {} \; | sort -h
> ```

### Desktop App

```bash
# Start Electron app (builds contracts, UI, starts Vite renderer + napi IPC bridge)
./dev/docker/run ./dev/dev d
```

### Web App

```bash
# Start web dev (Bacon auto-restarts Rust server on source changes, Vite serves frontend with HMR)
./dev/docker/run ./dev/dev a
```

## Testing

```bash
# All Rust tests
./dev/docker/run ./dev/dev t

# Specific tests by regex
./dev/docker/run ./dev/dev t gtr

# All TypeScript tests
./dev/docker/run ./dev/dev jt

# Full quality check (Rust lint + format + test)
./dev/docker/run ./dev/dev q

# Smoke tests
./dev/docker/run ./dev/dev br treetime && ./dev/docker/run ./dev/run-smoke-tests .out/treetime
```

## Linting and Formatting

```bash
# Rust
./dev/docker/run ./dev/dev l     # lint
./dev/docker/run ./dev/dev lf    # lint with auto-fix
./dev/docker/run ./dev/dev f     # format

# TypeScript (all packages)
./dev/docker/run ./dev/dev jl    # lint
./dev/docker/run ./dev/dev jlf   # lint with auto-fix
./dev/docker/run ./dev/dev jf    # format
```

## Monorepo Tooling

The TypeScript packages use [bun](https://bun.sh/) as the package manager and [turbo](https://turborepo.dev/) as the monorepo task runner.

Turbo handles the dependency chain automatically: starting the desktop app (`ds`) builds contracts and UI first. No manual `npm install` or intermediate build steps needed.

```bash
# Install all JS dependencies
./dev/docker/run ./dev/dev ji

# Upgrade all JS dependencies to latest versions
./dev/docker/run ./dev/dev ju
```

Configuration:

- `package.json` - bun workspace definition
- `turbo.json` - task dependencies and caching
- `bunfig.toml` - exact version pinning

## Performance Assessment

### Macro-benchmarks

```bash
cargo -q build --release --bin=treetime
export v='mpox/clade-ii/500'
hyperfine --warmup 1 --show-output "./target/release/treetime ancestral --method-anc=parsimony --outdir=tmp/$v --tree=data/$v/tree.nwk data/$v/aln.fasta.xz"
```

> 💡 Make sure you are benchmarking optimized (release) builds.

### Profiling

```bash
export v='mpox/clade-ii/500'
./dev/profile treetime ancestral --method-anc=parsimony --outdir="tmp/$v" --tree="data/$v/tree.nwk" "data/$v/aln.fasta.xz" -j1
```

> Requires additional configuration! Read comments inside the script first.

## Maintenance

### Upgrading Rust

The Rust version is defined in `rust-toolchain.toml`. When using `cargo`, the version defined in this file is installed automatically.

### Upgrading Rust Dependencies

```bash
./dev/docker/run ./dev/dev U
```

Note that dependency upgrades can cause breakage. The upgraded dependencies need to be reviewed, and unit tests, smoke tests and manual sanity checks may need to be performed.

### Upgrading JavaScript Dependencies

```bash
./dev/docker/run ./dev/dev ju
```

### Versioning

The workspace version in `Cargo.toml` tracks the base release version (e.g. `1.0.0`). Nightly builds append a prerelease suffix: `1.0.0-nightly.20260714T043012Z+a1b2c3d`.

### Releases

#### Stable releases

Not yet configured. The `dev/publish-github` script and the commented-out `publish-to-github-releases` job in `.github/workflows/cli.yml` are prepared for this.

#### Nightly releases

Automated pre-release builds are published in [neherlab/treetime](https://github.com/neherlab/treetime/releases) for early testing. Binaries for all 7 cross-compilation targets are attached to each release.

**Schedule**: daily at 04:00 UTC via `.github/workflows/schedule-nightly.yml` on the default branch. The dispatcher calls `.github/workflows/nightly.yml` on `rust` and skips if no new Rust commits exist since the last nightly.

**Manual trigger**:

```bash
./dev/trigger-nightly
```

This dispatches `schedule-nightly.yml` via `gh workflow run`. Requires `GH_TOKEN` or `gh auth login`.

**Version format**: `<cargo-version>-nightly.<YYYYMMDD>T<HHMMSS>Z+<short-sha>` (e.g. `1.0.0-nightly.20260714T043012Z+a1b2c3d`). Always marked as prerelease.

**How it works**:

1. `schedule-nightly.yml` on `master` calls `nightly.yml` on `rust`
2. `nightly.yml` resolves the current `rust` commit and calls `cli.yml` for the full build and validation matrix
3. `dev/publish-nightly` creates a prerelease in `neherlab/treetime` with all binaries attached

**Authentication**: the dispatcher grants `contents: write` to the repository's standard `GITHUB_TOKEN`. Nightly builds do not publish TreeTime Docker images; the reused CLI workflow retains its existing Docker builder-image cache.

### Continuous Integration (CI), Packaging, and Distribution

CI runs on every push to `rust` and on pull requests via `.github/workflows/cli.yml`:

- Cross-compilation for 7 targets (Linux gnu/musl, macOS, Windows)
- Unit tests
- Clippy lints
- Compatibility tests on native macOS, Windows, and Linux runners
- CLI documentation freshness check

The Rust nightly workflow (`.github/workflows/nightly.yml`) reuses `cli.yml` via `workflow_call` and adds a publish step.
