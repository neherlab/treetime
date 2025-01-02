# Contributing to Treetime

This guide describes how to set up a development environment, build Treetime, contribute to the codebase, and maintain the project.

> ⚠️ This is a document for developers, maintainers of the project as well as for the most curious users. It assumes basic familiarity with Treetime program from a user perspective. If you are not familiar with Treetime, read [User Documentation](TODO) first.

## Setup Development Environment

This guide assumes you are using Ubuntu 24.04, but the instructions will likely be similar for other Linux and Unix-like operating systems.

Treetime is written in Rust, and the standard `rustup` & `cargo` workflow is used.

## Install dependencies.

```bash
# Install required dependencies
# These commands are specific for Ubuntu Linux and may work on other Debian-based Linux distros.
# Refer to your operating system's documentation to find out how to install these dependencies.
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

# Add Rust tools to the $PATH. This line can be added to your .bashrc or .zshrc to adjust the $PATH automatically when a new terminal session is opened.
export PATH="$PATH:$HOME/.cargo/bin"
```

## Obtain Source Code

Treetime is an open-source project licensed under the MIT license, and its source code is available on GitHub.

To obtain the source code, use `git` to clone the GitHub repository:

```bash
git clone https://github.com/neherlab/treetime
```

If you are a team member, use the SSH URL to push securely and without a password prompt (More details: https://docs.github.com/en/authentication/connecting-to-github-with-ssh):

```bash
git clone git@github.com:neherlab/treetime.git
```

If you are not a team member but want to contribute, make a [fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo) and clone your forked repository instead. You can then submit changes to Treetime as a [Pull Request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request). Treetime maintainers will then review your changes and consider merging them into the main project.

## Build and Run

### Build and run the executables

```bash
# Go to the cloned directory
cd treetime

# (optional) checkout a branch (different from the default)
git checkout <branch_name>

# Build and run in debug mode (convenient for development, fast to build, slow to run, has more information in stack traces and when running under a debugger)
cargo run --bin=treetime

# Run with additional arguments passed to the executable
cargo run --bin=treetime -- --help

# Run Treetime in release mode (slow to build, fast to run, very little information in stack traces and during debugging)
cargo run --release --bin=treetime

# Alternatively, build and run separately. The compiled binaries will be in the `target/` directory by default.
# Step 1: build
cargo build --release --bin=treetime
# Step 2: run
./target/release/treetime

# (optional) To enable PNG image output, enable the "png" cargo feature
cargo run --features=png --bin=treetime
```

Quick examples of usage of different subcommands:

```bash
v="flu/h3n2/20"

# Run ancestral parsimony reconstruction on './data/flu/h3n2/20' 
cargo -q run --bin=treetime -- ancestral --method-anc=parsimony --outdir="tmp/ancestral-parsimony/$v" --tree="data/$v/tree.nwk" "data/$v/aln.fasta.xz"

# Run ancestral marginal sparse reconstruction on './data/flu/h3n2/20'
cargo -q run --bin=treetime -- ancestral --method-anc=marginal --model=jc69 --outdir="tmp/ancestral-marginal-sparse/$v" --tree="data/$v/tree.nwk" "data/$v/aln.fasta.xz"

# Run ancestral marginal dense reconstruction on './data/flu/h3n2/20'
cargo -q run --bin=treetime -- ancestral --method-anc=marginal --dense --model=jc69 --outdir="tmp/ancestral-marginal-sparse/$v" --tree="data/$v/tree.nwk" "data/$v/aln.fasta.xz"

# Run clock estimation on './data/flu/h3n2/20'
cargo -q run --bin=treetime -- clock --tree="data/$v/tree.nwk" --dates="data/$v/metadata.tsv" --outdir="tmp/clock/$v"
```

> 💡 Set variable `v` (virus) to a different path in the `data/` directory to try the same commands with other example inputs. List names of all examples with:
> 
> ```
> find data -type f -name "tree.nwk" -exec dirname {} \; | sort -h
> ```
> 
> Note that some inputs are very large and may take very long time to run, especially in debug mode.

## Testing

### Unit tests

Tests are run using [nextest](https://nexte.st/). This can be installed with:
```bash
cargo install nextest
```

Then run the tests with:
```bash
cargo nextest run
```

Add the `--no-fail-fast` flag to continue running tests even if there are failures.

You can run a subset of tests by providing a regex matching the full test name. For example:

```bash
cargo nextest run gtr
```

See also: [nextest running tests](https://nexte.st/docs/running/)


### Smoke tests

Smoke tests run different Treetime commands on many example inputs in bulk. This serves as a screening for bugs and other regressions. Note that the script does not rebuild Treetime - you need to (re-)build the executable yourself and then pass the path to it to the script:

```bash
cargo run --release --bin=treetime
./dev/run-smoke-tests ./target/release/treetime
```

Refer to the comments in the script for details. Set variable `verbose = True` in the script to print the commands being run.


## Linting (Static Analysis)

Rust code is linted by running [Clippy](https://github.com/rust-lang/rust-clippy):

```bash
cargo clippy
```

Apply automatic fixes

> ⚠️ Destructive operation! Save your code changes first!

```bash
cargo clippy --fix
```

Clippy is configured in `clippy.toml` and `.cargo/config.toml`.

## Formatting (Code Style)

Code formatting is handled by [rustfmt](https://rust-lang.github.io/rustfmt):

> ⚠️ Destructive operation! Save your code changes first!

```bash
cargo fmt --all
```

Rustfmt is configured in `rustfmt.toml`.

## Performance assessment

### Macro-benchmarks

In order to measure running time of a CLI program end-to-end, you could use [`hyperfine`](https://github.com/sharkdp/hyperfine).

Example:

```bash
cargo -q build --release --bin=treetime
export v='mpox/clade-ii/500'
hyperfine --warmup 1 --show-output "./target/release/treetime ancestral --method-anc=parsimony --outdir=tmp/$v --tree=data/$v/tree.nwk data/$v/aln.fasta.xz"
```

Example output:
```
  Time (mean ± σ):     655.4 ms ±  25.7 ms    [User: 2258.9 ms, System: 210.1 ms]
  Range (min … max):   635.3 ms … 717.6 ms    10 runs
```

> 💡 Make sure you are benchmarking optimized (release) build of a program.

> 💡 To obtain representative and reproducible benchmarking results, make sure you are benchmarking on a "calm" system with minimum number of concurrent processes and with no shortage of memory and disk space etc.

### Profiling

Sampling (statistical) performance profiling can be performed using helper script `./dev/profile`:

> Requires additional configuration! Read comments inside the script first!

```bash
export v='mpox/clade-ii/500'
./dev/profile treetime ancestral --method-anc=parsimony --outdir="tmp/$v" --tree="data/$v/tree.nwk" "data/$v/aln.fasta.xz" -j1
```

## Maintenance

### Upgrading Rust

The Rust version is defined in `rust-toolchain.toml`. When using `cargo`, the version defined in this file is installed automatically.

### Upgrading Dependencies

Dependencies for subprojects are defined in `packages/**/Cargo.toml` and in `Cargo.lock`. They are periodically upgraded by a dedicated maintainer, manually using `cargo-upgrade` from the [cargo-edit](https://github.com/killercup/cargo-edit) package.

```bash
export CARGO_NET_GIT_FETCH_WITH_CLI=true
export CARGO_REGISTRIES_CRATES_IO_PROTOCOL='git'
export CARGO_REGISTRIES_CRATES_IO_PROTOCOL='git cargo fetch'
export CARGO_REGISTRIES_CRATES_IO_PROTOCOL='git cargo upgrade'
cargo upgrade --pinned --incompatible --verbose
```

Note that dependency upgrade can cause serious breakage - both compile-time and run-time. The upgraded dependencies, including sub-dependencies, need to be reviewer, unit tests, smoke tests and some manual sanity checks may need to be performed.

### Versioning

TODO

### Releases

TODO

### Continuous Integration (CI), Packaging, and Distribution

TODO
