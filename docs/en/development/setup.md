---
title: Setup
description: Build TreeTime in the build container or on the host, and the optional settings of a machine.
---

TreeTime builds in two ways: in the build container (recommended), or directly on the host. Both run the same `just` recipes with the same tool versions, which `.config/mise.toml` pins and `.config/mise.lock` locks by URL and checksum.

## Build container (recommended)

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

## Host

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

## Machine settings

Optional settings go into the gitignored `.env` in the checkout; `.env.example` lists them:

- `KACHE_STORE`: directory of the [kache](https://github.com/kunobi-ninja/kache) compiler cache stores. Every build and clippy pass then compiles through kache, which shares compiled crates across the worktrees of this project. Workspace crates of incremental builds (`dev`, tests, `release`, clippy) keep their incremental state, so an edit rebuilds as fast as without kache; kache serves their dependencies and every build without incremental state (CI, `dist`, `profiling`, `bench`, cross builds). Dylint, hawk, and coverage compile without kache. The directory holds one store per kache version, environment (`host` or `docker`), and pass (`build`, `clippy`, or `cross-<target>`)
- `KACHE_MAX_SIZE`: size limit of each kache store, 100 GiB when unset
- `TREETIME_DOCKER_CPUS`: CPU budget of the containers that the dev scripts start, 8 when unset and at most the CPUs of the Docker daemon. A container of `dev/docker/run` uses at most this many CPUs; `dev/cross/all` and the Linux distribution tests divide the budget between the containers they run at the same time. Cargo in a container runs one job per whole CPU of its budget unless `CARGO_BUILD_JOBS` is set, so the compiler work stays within the budget also where a host daemon such as ananicy moves compilers out of the container cgroup
- `SMOKE_STORE`: absolute directory of a smoke snapshot store shared by all checkouts of this project on the machine, so each baseline is built once (see [Snapshot store](/development/testing#snapshot-store)). Unset: each checkout keeps its own `snapshots/`
- `SMOKE_MIN_FREE_GIB`: free disk space that `dev/smoke` keeps on the file systems of the checkout and the snapshot store, 20 GiB when unset (see [Disk limits](/development/testing#disk-limits))
