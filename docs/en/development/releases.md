---
title: Cross-compilation and releases
description: Cross builds, nightly releases, Docker images, and the deploys of the web app and the docs.
---

`just cross` (host) builds the shipped CLI (`dist` profile) for every target in `dev/cross/targets` in its cross image, in parallel, into `.out/`. The targets divide the CPU budget `TREETIME_DOCKER_CPUS` between them (see [Machine settings](/development/setup#machine-settings)). `just cross --target=aarch64-apple-darwin` builds one target. Without `just` on the host, run `./dev/cross/all treetime`.

`just cross-napi` (host) builds the Node addon of the desktop app the same way, with the same profile and CPU flags, for every target in `dev/cross/targets-desktop`, the targets Electron supports (no musl), into `.out/app-napi-<target>.node` (`./dev/cross/all --lib app-napi` without `just`). The Windows addon is built with llvm-mingw, like the CLI, and loads into Electron, which is built with MSVC: Node-API reaches the addon through the functions that `napi-sys` looks up in the running executable, not through an import library.

Each target family has its own cross image. mise installs every toolchain component and library archive of the images from `.config/mise.image.toml` at the URL and checksum of `.config/mise.image.lock` (`dev/docker/files/install-image-tools`), and the libraries link through pkg-config, as in the native image (see [Build container](/development/setup#build-container-recommended)):

- Linux (`dev/docker/cross-linux.dockerfile`): clang and lld of LLVM 22 ([binarylandia/build_llvm](https://github.com/binarylandia/build_llvm)) on the sysroot, crt files, libgcc, and libstdc++ of a GCC 16 toolchain of [binarylandia/build_crosstool-ng](https://github.com/binarylandia/build_crosstool-ng): glibc 2.17 for the gnu targets, musl 1.2.6 for the musl targets. The musl binaries are static
- Windows (`dev/docker/cross-windows.dockerfile`): llvm-mingw of the same LLVM release, with the static mingw-w64 sysroot on the Universal C Runtime (UCRT), compiler-rt, and libunwind, for `x86_64-pc-windows-gnullvm` and `aarch64-pc-windows-gnullvm`, and the OpenBLAS, liblzma, and libzstd archives built with it. The binaries run on Windows 10 and later
- Windows 7 (`dev/docker/cross-windows7.dockerfile`, target `x86_64-win7-windows-gnu`, not in `dev/cross/targets`): the GCC 16 MinGW toolchain of build_crosstool-ng on the msvcrt C runtime, which Windows 7 has without updates, the msvcrt OpenBLAS of build_openblas, and the liblzma and libzstd archives of staticbourg built with this GCC toolchain. Rust ships no standard library for this tier 3 target, so `dev/cross/build` compiles it from the `rust-src` component with `-Zbuild-std` (and `RUSTC_BOOTSTRAP=1`, which the stable toolchain needs for that flag). The CLI only: Electron does not run on Windows 7. Build it with `just cross --target=x86_64-win7-windows-gnu`
- WebAssembly (`dev/docker/cross-wasm.dockerfile`, target `wasm32-unknown-emscripten`, not in `dev/cross/targets`): Emscripten, the Rust toolchain with the nightly of the threaded module, the WebAssembly archives of OpenBLAS, liblzma, and libzstd, and the tools of `just build-web2`, without the native toolchains, the lint toolchains, and the other development tools of the native image. The GCC of Debian compiles the build scripts. The web deploy builds in it, and so can you: `CROSS_COMPILE=wasm32-unknown-emscripten ./dev/docker/run just build-web2`. The development image keeps Emscripten, for the WebAssembly lints and tests of `just check`
- macOS (`dev/docker/cross-darwin.dockerfile`): the llvm flavor of osxcross ([binarylandia/build_osxcross](https://github.com/binarylandia/build_osxcross)) with the macOS 11.1 SDK, which runs clang and `ld64.lld` of the LLVM release, and compiler-rt. The deployment targets are macOS 10.12 (x86_64) and 11.0 (arm64)

`dev/cross/check` fails when a binary or addon needs a library that is not part of the operating system: anything beyond glibc 2.17 on the Linux gnu targets (the addon may also need `libgcc_s.so.1`, as Electron does), any library on the musl targets, a DLL other than the Windows system DLLs (`KERNEL32.dll`, `api-ms-win-crt-*.dll`, and others) on Windows, a DLL other than the Windows 7 system DLLs and `msvcrt.dll`, or a function that Windows 8 or later added, on Windows 7, and a library outside `/usr/lib` and `/System/Library` on macOS. It also fails when a macOS binary names another minimum macOS version than its deployment target. `dev/cross/all --run` runs each Linux binary under QEMU and each x86_64 Windows binary under Wine; macOS and aarch64 Windows binaries are only built and checked there. For each Linux binary it also runs `dev/cross/check-cpus`, which runs the binary under QEMU on the oldest CPU of its baseline (Nehalem, Cortex-A53) with every compression format, and on x86_64 below the baseline, where the binary must stop with its CPU message. The nightly and a manual run of `CI` run `dev/cross/check` after each build, and `.github/workflows/cli-compat.yml` runs `dev/cross/check-cpus` in its `cpu-baseline` jobs, the Linux binaries on distribution images for x86_64 and arm64 (`compat-linux`), the Windows binaries on `windows-2022` and `windows-11-arm` (`compat-windows`), and the macOS binaries on Intel and Apple Silicon macOS runners (`compat-macos`).

The release binaries run on these CPUs; on an older x86_64 CPU they stop with a message that names the missing level. The prebuilt archives they link follow the same levels: build_openblas and staticbourg compile them with the C CPU flags of `export_dist_flags`, which the image build checks in the build record of every archive, and OpenBLAS also selects its kernels at run time.

- x86_64 (every operating system, Windows 7 included): x86-64-v2 (SSE4.2, POPCNT, CMPXCHG16B; Intel Nehalem, AMD Bulldozer, or newer)
- Linux aarch64: Armv8.0, such as the Cortex-A53 and Cortex-A72 of the Raspberry Pi 3 and 4. Atomic operations go through outline atomics, which use LSE instructions when the CPU has them
- Windows aarch64: the LSE atomics of Armv8.1, which every device that runs Windows 11 on Arm has
- macOS aarch64: Apple M1 or newer

## Nightly releases

Prerelease builds are published in [neherlab/treetime-nightly](https://github.com/neherlab/treetime-nightly/releases). The nightly is the only workflow that builds the CLI on `rust`, and it deploys every binary that builds: checks and tests report problems but never hold back a deployment. No workflow builds the desktop app; `just package-desktop` packages it locally (see [Desktop app](/development/desktop-app)).

1. `.github/workflows/nightly.yml` on `rust` runs on every push to `rust` and releases its tip. Runs do not overlap: a push during a running nightly waits, and only the newest waiting run starts, so one release can cover several commits. `.github/workflows/schedule-nightly.yml` on `master` also calls `nightly.yml` daily at 04:00 UTC, which skips when the latest nightly already has the tip of `rust`, so it releases only a commit whose push run failed. The scheduler passes the inputs `force` and `only`; `nightly.yml` accepts both and ignores `only`, so the scheduler on `master` works unchanged
2. `nightly.yml` builds the cross-compilation matrix of `.github/workflows/cli-build.yml`. A failed target does not block the others. Each binary is uploaded before `dev/cross/check` checks its linkage, so a linkage problem is reported and the binary still ships
3. `nightly.yml` writes `SHA256SUMS` over every asset of the release and attests these checksums with `actions/attest`, which signs with the identity of the workflow run. The attestation needs the `id-token: write` and `attestations: write` permissions. A called workflow can only lower the permissions of its caller, so the `nightly` job of `schedule-nightly.yml` grants them too
4. `dev/publish-nightly` creates the prerelease with the built binaries, `examples.zip`, which holds the files that `data/index.json` lists and the index itself, and `SHA256SUMS`. Its notes list the changelog since the previous nightly, followed by a downloads table: the CLI binaries, one row per platform, the `docker pull` command of the nightly image, the other assets, and the command that verifies a download: `gh attestation verify <file> --repo neherlab/treetime`
5. In parallel with the release, the `docker` job publishes the Docker images `neherlab/treetime` with `dev/docker/prod-images --push --run`, from the Linux binaries of the same run (see [Docker images](#docker-images))
6. In parallel with the release, `.github/workflows/cli-compat.yml` runs the same binaries on Linux distributions, on the oldest CPU of each Linux baseline, on macOS, and on Windows
7. After the release, `.github/workflows/nextstrain-runtimes.yml` starts the CI of [nextstrain/docker-base](https://github.com/nextstrain/docker-base) and [nextstrain/conda-base](https://github.com/nextstrain/conda-base), whose builds download the latest release as `treetime1` and publish the Nextstrain runtimes (the image `nextstrain/base` and the Conda package `nextstrain-base`). It skips a repository whose build is queued but not yet started, because that build downloads the new release when it starts, and it never waits for the builds. It uses `GH_TOKEN_NEXTSTRAIN_BOT_WORKFLOW_DISPATCH_NEHERLAB`, a classic token of the Nextstrain bot account with the `public_repo` scope, stored as an organization secret of neherlab that only this repository can read

`./dev/trigger-nightly` starts a nightly by hand; `--force` builds also when the latest nightly already has the tip of `rust`. The version format is `<cargo-version>-nightly.<YYYYMMDD>T<HHMMSS>Z+<short-sha>`.

## Docker images

`dev/docker/prod-images --version-suffix <suffix> [--run] [--push]` (host only) builds the images `neherlab/treetime` for `linux/amd64` and `linux/arm64` from the four Linux CLI binaries of a nightly in `.out/` (`treetime-{x86_64,aarch64}-unknown-linux-{gnu,musl}`). The images copy the binaries and compile nothing:

- **Debian** (`dev/docker/prod-debian.dockerfile`): the glibc build on the newest Debian, with `bash`, `ca-certificates`, `curl`, `procps`, and `wget`, which workflow managers such as Nextflow need in task containers
- **Alpine** (`dev/docker/prod-alpine.dockerfile`): the static musl build on the newest Alpine
- **Scratch** (`dev/docker/prod-scratch.dockerfile`): only the static musl build

`nightly-debian`, `nightly-alpine`, and `nightly-scratch` follow the latest nightly, and `nightly` is the Debian image. The version of a nightly selects one nightly, with `_` in place of the `+` that a Docker tag cannot hold: `1.0.0-nightly.20261026T040000Z_abc1234` (Debian) and `1.0.0-nightly.20261026T040000Z_abc1234-alpine`. There is no `latest` tag until stable releases exist. Each image is built for every platform whose binary is in `.out/`, so a missing binary leaves out only its platform. With `--push --run`, the script pushes every image first, then loads each platform and runs `treetime --version` and an ancestral reconstruction of `data/zika/20` in every image. A missing binary, a failed build, or a failed test fails the script at the end, after every other image is pushed and tested. Usage:

```sh
docker run --rm --user="$(id -u):$(id -g)" --volume="$PWD:/data" --workdir=/data neherlab/treetime:nightly treetime ancestral --tree=tree.nwk --alignment=aln.fasta.xz --output-all=results
```

`dev/docker/prod-bases.json` pins the Debian and the Alpine base to the digest of its Docker Hub tag, so the images of one nightly and a re-run use the same base. `./dev/docker/run dev/docker/prod-bases-update` rewrites it from the release lines of endoflife.date and the tags of Docker Hub, with the seven-day rule of the project:

- the newest release line released at least seven days ago is used
- when the tag of the line was pushed less than seven days ago, the line takes the newest dated tag of the same line that is at least seven days old: Debian keeps every build as `<codename>-YYYYMMDD`, such as `trixie-20260918`; Alpine has no dated tags, so the Alpine line keeps the digest that the file already has
- a line without such a digest fails the run, with the date from which a re-run succeeds

Run it from time to time and commit the file, so the images get the security fixes of their bases. The repository page on Docker Hub is edited by hand. The `docker` job logs in with `NEHERLAB_BOT_DOCKERHUB_USERNAME` and `NEHERLAB_BOT_DOCKERHUB_TOKEN`, the account that also pushes the build images, which needs write access to `neherlab/treetime`.

## Web app deploys

`.github/workflows/web.yml` builds the web2 app (`just build-web2`) on every push in the WebAssembly image (`CROSS_COMPILE=wasm32-unknown-emscripten`, see the WebAssembly image above) and uploads it with `dev/publish-web` to the preview host, a CloudFront distribution in front of an S3 bucket in the AWS account `preview-deployments` of neherlab. `just build-web2` compiles the single-threaded and the threaded module at the same time. `dev/publish-web` uploads the build to the first site and copies it to the other sites inside the bucket. Each upload goes to several sites:

- `https://neherlab--treetime--commit-<short-sha>.previews.neherlab.click`: the commit
- `https://neherlab--treetime--branch-<branch>.previews.neherlab.click`: the latest upload of the branch
- `https://neherlab--treetime--pr-<number>.previews.neherlab.click`: the latest upload of each open pull request of the branch; the run adds or updates a comment with the link in the pull request
- `https://treetime.neherlab.click`: the latest push to `rust`

Opening or reopening a pull request also starts a run, so a pull request gets its site without a new push. Pull requests from forks build without an upload, because workflows of forks get no secrets. The upload reads the access key, the region, the bucket, and the distribution from the organization secrets `PREVIEWS_NEHERLAB_CLICK_AWS_ACCESS_KEY_ID`, `PREVIEWS_NEHERLAB_CLICK_AWS_SECRET_ACCESS_KEY`, `PREVIEWS_NEHERLAB_CLICK_AWS_DEFAULT_REGION`, `PREVIEWS_NEHERLAB_CLICK_AWS_S3_BUCKET`, and `PREVIEWS_NEHERLAB_CLICK_AWS_CLOUDFRONT_DISTRIBUTION_ID`. `dev/publish-web --help` describes the stored files, their caching, and what the edge functions of the distribution do; the edge functions themselves are not in this repository.

## Docs deploys

`.github/workflows/docs.yml` builds the docs site (`just build-docs`) on every push and uploads it with the same composite action as the web app (`.github/actions/publish-preview`). Each docs build is a version folder of its site, today only `dev`, the tip of `rust`:

- `https://neherlab--treetime-docs--commit-<short-sha>.previews.neherlab.click/dev`, `...--branch-<branch>...`, and `...--pr-<number>...` (with a pull request comment marked `<!-- docs-preview -->`)
- `https://docs.treetime.neherlab.click/dev`: the latest push to `rust`, from the bucket folder `treetime-docs`

Each upload also stores the site root files of `dist-root/`: `versions.json`, the list of versions, and an `index` page that leads to the preferred version. The version select and the banner of every build read `/versions.json` at run time, so an older build learns about a newer version without a rebuild. `packages/app-docs/versions.json` is the source of the list; `TREETIME_DOCS_VERSION` names the version of a build (default `dev`), and the build fails when the list has no entry for it. Adding a release version means building the docs at the release tag with `TREETIME_DOCS_VERSION=<version>`, uploading with `dev/publish-web --prefix <version>`, and adding the entry, marked `preferred`, to `versions.json`.
