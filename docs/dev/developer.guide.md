# Contributing to Treetime

This guide describe how to setup developer environment, how to build Treetime, contribute to the codebase and maintain
the project.


## Setup developer environment

This guide assumes Ubuntu 24.04 operating system, but will likely work similarly to any other Linux and Unix-like
machine.

Treetime is written in Rust. The usual `rustup` & `cargo` workflow can be used:

```bash
# Install required dependencies
# These particular commands are specific for Ubuntu Linux and will work on some other Debian-based Linux distros.
# Refer to documentation of your operating system to find how to install these dependencies.
sudo apt-get update
sudo apt-get install bash clang curl gcc gfortran git make pkg-config protobuf-compiler libopenblas-dev

# (optional) if you want to enable png image output ("png" cargo feature, see below), then add
sudo apt-get install libfontconfig1-dev

# Install Rustup, the Rust version manager (https://www.rust-lang.org/tools/install)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y

# Add Rust tools to the $PATH. You might add this line to your .bashrc or .zshrc so that the $PATH is adjusted automatically when a new terminal session is opened.
export PATH="$PATH:$HOME/.cargo/bin"

````

## Obtain source code

Treetime is an open-source project and its source code is available on GitHub under MIT license.

To obtain source code use `git` to clone GitHub repository:

```bash
git clone https://github.com/neherlab/treetime
```

If you are a team member, use SSH url to be able to push securely and without password prompt
(More details: https://docs.github.com/en/authentication/connecting-to-github-with-ssh)

```bash
git clone git@github.com:neherlab/treetime.git
```

If you are not a team member, but want to contribute, make a
[fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo), and clone your forked repository instead. You can
then submit changes to treetime as
a [Pull Request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request)
. Treetime maintainers will then review your changes and will consider merging them into the main project.

## Build and run

### Build and run executables

```bash
# Go to the cloned directory
cd treetime

# (optional) checkout a branch (different from default)
git checkout <branch_name>

# Build and run in debug mode (convenient for development, fast to build, slow to run, has more information in stack traces and when running under a debugger)
cargo run --bin=treetime

# Instead of `--bin=treetime` you can also run any other executable from `packages/treetime/src/bin/`. Just substitute its filename.
# This is a `cargo` convention everything in `src/bin/` that has a `main()` function in it becomes an executable. This way you can add more executables.
cargo run --bin=graph
cargo run --bin=sequence_algorithms

# Run Treetime in release mode (slow to build, fast to run, very little information in stack traces and during debugging)
cargo run --release --bin=treetime

# Alternatively, build and run separately. The compiled binaries will be in `target/` directory by default.
cargo run --release --bin=treetime
./target/release/treetime

# (optional) To allow png image output, then enable "png" cargo feature
cargo run --bin=treetime --feature=png
```

Note, on first build of a particular project, `cargo` will search for one of the possible toolchain config files and
will automatically install Rust version required by the project. This may cause first build to take longer than usual.

### Testing

Run all test with:

```bash
cargo test
```

Add `--no-fail-fast` flag to keep going even if there are failures.

A subset of tests can be ran by providing a regex matching full test name. For example

```bash
cargo test gtr
```

See also: [cargo-test](https://doc.rust-lang.org/cargo/commands/cargo-test.html)

### Linting (static analysis)

Rust code is linted by running [Clippy](https://github.com/rust-lang/rust-clippy):

```bash
cargo clippy
```

Clippy is configured in `clippy.toml` and `.cargo/config.toml`.

### Formatting (code style)

Code formatting is done using [rustfmt](https://rust-lang.github.io/rustfmt):

```bash
cargo fmt --all
```

Rustfmt is configured in `rustfmt.toml`.

## Maintenance

### Upgrading Rust

Rust version is defined in `rust-toolchain.toml`. When using `cargo`, the version defined in this file gets installed
automatically.

### Upgrading Rust

Dependencies for subprojects are defined in  `packages/**/Cargo.toml` and in `Cargo.lock`. They are periodically
upgraded by a dedicated maintaier, manually using `cargo-upgrade`
from [cargo-edit](https://github.com/killercup/cargo-edit) package.

```bash
cargo upgrade --workspace
```

### Versioning

TODO

### Releases

TODO

### Continuous integration (CI), packaging and distribution

TODO
