# TreeTime task runner: `just` lists every recipe by group.
#
# Recipes run the tools pinned in mise.toml, on the host after `mise install` or
# in the build container through `dev/docker/run just <recipe>`. Before a change
# is merged, `just check-all` must pass.
#
# Naming: a leaf recipe runs one tool in one mode. The suffix `-rs` or `-ts`
# names the toolchain of a leaf, Rust or TypeScript (Bun); tools that serve one
# toolchain only, such as dylint or knip, keep their own name. A recipe without
# a language suffix combines the leaves of both toolchains and calls no tool
# itself. The suffix `-all` adds the slow tools to the fast set of the same name:
# `fix` and `fix-all`, `lint` and `lint-all`, `test` and `test-all`, `check` and
# `check-all`.

set minimum-version := "1.58.0"
set default-list
set dotenv-load
set positional-arguments
set shell := ["mise", "exec", "--", "bash", "-euo", "pipefail", "-c"]
set script-interpreter := ["mise", "exec", "--", "bash", "-euo", "pipefail"]

export CARGO_TERM_QUIET := "true"
export NEXTEST_NO_TESTS := "fail"

# Build output of the checkout. dev/docker/run points TREETIME_BUILD_DIR at
# .build/container, so host and container artifacts, which link against
# different system libraries, never mix. Builds and tests share one cargo target
# directory; clippy, dylint, and hawk each have their own, so a lint never waits
# on the build lock.
build_dir := env("TREETIME_BUILD_DIR", justfile_directory() / ".build/host")
export CARGO_TARGET_DIR := build_dir / "cargo"
lint_target_dir := build_dir / "lint"
dylint_target_dir := build_dir / "dylint"
hawk_target_dir := build_dir / "hawk"

# The kache compiler cache is optional. KACHE_STORE, in .env or the environment,
# names the store directory and switches builds and clippy to kache, each pass
# with its own store. Dylint and hawk always compile without it, because a cache
# hit skips the lint passes.
kache_store := env("KACHE_STORE", "")
export RUSTC_WRAPPER := if kache_store != "" { "kache" } else { env("RUSTC_WRAPPER", "") }
export KACHE_CACHE_DIR := if kache_store != "" { kache_store / "build" } else { env("KACHE_CACHE_DIR", "") }
lint_env := "CARGO_TARGET_DIR=" + quote(lint_target_dir) + if kache_store != "" { " KACHE_CACHE_DIR=" + quote(kache_store / "lint") } else { "" }
uncached_env := "RUSTC_WRAPPER= CARGO_INCREMENTAL=0"

# The dylint driver, shared by the check, fix, and baseline recipes. The
# pub_unused_in_workspace lint writes its records into the pub-unused directory,
# and pub-unused-report reads them after the check pass. Mordant lists the
# crates over the committed baseline in over-baseline.txt.
pub_unused_env := "TREETIME_LINTS_PUB_UNUSED_DIR=" + quote(dylint_target_dir / "pub-unused")
dylint_cmd := uncached_env + " CARGO_TARGET_DIR=" + quote(dylint_target_dir) + " " + pub_unused_env + " cargo dylint --quiet --all"
dylint_cargo_args := "--quiet --locked --workspace --all-targets"
dylint_over_baseline := dylint_target_dir / "mordant/over-baseline.txt"

# Library crates whose public API is an external boundary, skipped by both
# unused-public-code checks (hawk and pub-unused-report): the shared utility and
# file-format libraries publish a complete API for reuse, and the Node addon's
# `#[napi]` surface is consumed by JavaScript. The test-infrastructure crates
# serve only test targets, which neither check counts as users.
public_api_crates := "treetime_utils util_newick util_phyloxml util_augur_node_data_json util_usher_mat app_napi treetime_validation treetime_analytical"

# Cargo with the seven-day minimum publish age of new dependency releases. The
# setting is unstable in the pinned Rust (stable from 1.100), so the dependency
# recipes enable it for their resolution only; RUSTC_BOOTSTRAP=1 allows unstable
# Cargo features on a stable toolchain.
cargo_min_age := "RUSTC_BOOTSTRAP=1 cargo -Zmin-publish-age --config 'registry.global-min-publish-age=\"7 days\"'"
hawk_toolchain := trim(read("dev/docker/files/hawk-toolchain"))

# Groups of the full gate, one CI job each (`just check-group <group>`). Hawk
# needs a second full compilation and runs locally only.
checks_format := "fmt-check-rs fmt-check-ts fmt-check-other lint-shell lint-docker lint-workflows deny shear fixtures-check"
checks_clippy := "lint-rs"
checks_dylint := "dylint test-dylint"
checks_tests := "test-rs"
checks_generated := "generated-check"
checks_typescript := "typecheck lint-ts knip test-ts"
checks_hawk := "hawk"
check_fast := "fmt-check-rs fmt-check-ts fmt-check-other lint-rs lint-ts typecheck"
check_full := checks_format + " " + checks_clippy + " " + checks_dylint + " " + checks_tests + " " + checks_generated + " " + checks_typescript + " " + checks_hawk

alias b := build
alias br := build-release
alias bp := build-profiling
alias r := run
alias rr := run-release
alias t := test-rs
alias tu := test-unit-rs
alias ti := test-integration-rs
alias sm := smoke
alias l := lint-rs
alias lf := lint-fix
alias f := fmt
alias fc := fmt-check

# Fast checks: formatting, clippy, TypeScript types and lints, in parallel
[group("check")]
check: _js
    TREETIME_JS_READY=1 dev/run-checks {{ check_fast }}

# Every check, in parallel; must pass before a change is merged
[group("check")]
check-all: _js
    TREETIME_JS_READY=1 dev/run-checks {{ check_full }}

# One group of check-all, serially with streamed output, as its CI job runs it: just check-group <format|clippy|dylint|tests|generated|typescript|hawk>
[group("check")]
check-group group: _js
    recipes="$(just --evaluate "checks_$1")"; TREETIME_JS_READY=1 dev/run-checks --serial ${recipes}

# Apply the fast automatic lint fixes (clippy, oxlint), then format; stage your changes first
[group("check")]
fix: lint-fix fmt

# Apply every automatic lint fix, dylint included, then format; stage your changes first
[group("check")]
fix-all: lint-fix dylint-fix fmt

# Install the pinned tools and lint toolchains (on the host) and the JavaScript dependencies
[group("setup")]
setup:
    if [[ -z "${TREETIME_CONTAINER:-}" ]]; then mise install; for dir in dev/lints/*/; do (cd "${dir}" && rustup toolchain install); done; rustup toolchain install {{ hawk_toolchain }} --profile minimal --component rustc-dev,llvm-tools-preview,rust-src; fi
    just _js

# Build the workspace (dev profile)
[group("build")]
build *args:
    cargo build --locked "$@"

# Build a binary (release profile) and copy it to .out/: just build-release treetime
[group("build")]
build-release bin="treetime" *args:
    cargo build --locked --release --bin {{ quote(bin) }} "${@:2}"
    mkdir -p .out && cp {{ quote(CARGO_TARGET_DIR / "release" / bin) }} .out/

# Build a binary (profiling profile: release with debug symbols)
[group("build")]
build-profiling bin="treetime" *args:
    cargo build --locked --profile=profiling --bin {{ quote(bin) }} "${@:2}"

# Cross-compile release binaries in the cross images (host only, needs Docker): just cross [--target=<triple>]
[group("build")]
cross *args:
    dev/cross/all "$@" treetime

# Run a binary (dev profile): just run treetime ancestral --help
[group("run")]
run bin *args:
    args=("${@:2}"); [[ "${args[0]:-}" != "--" ]] || args=("${args[@]:1}"); cargo run --locked --bin {{ quote(bin) }} -- ${args[@]+"${args[@]}"}

# Run a binary (release profile): just run-release treetime ancestral --help
[group("run")]
run-release bin *args:
    args=("${@:2}"); [[ "${args[0]:-}" != "--" ]] || args=("${args[@]:1}"); cargo run --locked --release --bin {{ quote(bin) }} -- ${args[@]+"${args[@]}"}

# Run an example (release profile): just example validation_test
[group("run")]
example name *args:
    cargo run --locked --release --example {{ quote(name) }} -- "${@:2}"

# Rust and TypeScript tests
[group("test")]
test: test-rs test-ts

# Rust and TypeScript tests, and the tests of the custom dylint libraries
[group("test")]
test-all: test test-dylint

# Rust tests (nextest); arguments are nextest filters and options
[group("test")]
test-rs *args:
    cargo nextest run --locked --workspace "$@"

# Rust unit tests only
[group("test")]
test-unit-rs *args:
    cargo nextest run --locked --workspace --lib "$@"

# Rust integration tests only
[group("test")]
test-integration-rs *args:
    cargo nextest run --locked --workspace --test '*' "$@"

# List Rust tests without running them
[group("test")]
test-list-rs *args:
    cargo nextest list --locked --workspace "$@"

# TypeScript tests (vitest) and the custom oxlint rule tests
[group("test")]
test-ts: _js
    bun run --silent test
    node --test "dev/lints/oxlint/__tests__/test_*.ts" "dev/lints/oxlint-anti-slop/**/*.test.ts"

# Tests of the custom dylint libraries and the pub-unused-report tool
[group("test")]
test-dylint:
    cd dev/lints/dylint-custom && {{ uncached_env }} cargo test --quiet --release --locked --lib --bins --target-dir {{ quote(dylint_target_dir / "report") }}
    cd dev/lints/dylint-trailofbits && {{ uncached_env }} cargo test --quiet --locked --workspace --lib --target-dir {{ quote(dylint_target_dir / "trailofbits-test") }}

# Smoke-test the fast cases (datasets of at most 100 sequences) against the rust branch (host only, needs Docker): just smoke [--only REGEX] [--against REF]
[group("test")]
smoke *args:
    dev/smoke "$@"

# Smoke-test all cases, every dataset and the slow extra commands, against the rust branch (host only, needs Docker)
[group("test")]
smoke-all *args:
    dev/smoke --tier full "$@"

# Smoke-test the CLI without a baseline: crash, timeout and output checks only (host only, needs Docker)
[group("test")]
smoke-run *args:
    dev/smoke --no-compare "$@"

# Smoke-test again the cases that did not pass in the last run (host only, needs Docker)
[group("test")]
smoke-failed *args:
    dev/smoke --rerun-failed "$@"

# Delete the smoke snapshots of dirty working trees other than the current one
[group("test")]
smoke-prune:
    dev/smoke --prune

# Rust and TypeScript coverage reports, never a gate
[group("report")]
coverage: coverage-rs coverage-ts

# Rust coverage report (cargo-llvm-cov), never a gate
[group("report")]
coverage-rs:
    {{ uncached_env }} cargo llvm-cov nextest --locked --workspace --html --output-dir tmp/coverage/rs
    printf 'Rust coverage: tmp/coverage/rs/html/index.html\n'

# TypeScript coverage report (vitest), never a gate
[group("report")]
coverage-ts: _js
    bun run --silent coverage
    printf 'TypeScript coverage: tmp/coverage/ts/index.html\n'

# Mutation testing of the code changed since the fork point (or since `rust`), never a gate
[group("report")]
mutants *args:
    dev/mutants "$@"

# Copy-paste duplication across Rust and TypeScript (jscpd), never a gate
[group("report")]
duplication *args:
    jscpd --config .jscpd.json "$@"

# Inventory of lint suppressions and other review-sensitive settings
[group("report")]
review-suppressions:
    dev/review-suppressions

# Fast lints: clippy and oxlint
[group("lint")]
lint: lint-rs lint-ts

# Every lint: the fast lints, the custom lint libraries, unused code and dependencies, and the shell, Dockerfile, and workflow lints
[group("lint")]
lint-all: lint dylint hawk deny shear knip lint-shell lint-docker lint-workflows

# Apply the fast automatic lint fixes: clippy, then oxlint; stage your changes first
[group("lint")]
lint-fix: lint-fix-rs lint-fix-ts

# Clippy over all targets, denying warnings
[group("lint")]
lint-rs *args:
    {{ lint_env }} cargo clippy --locked --workspace --all-targets "$@" -- --deny warnings

# Apply clippy's machine-applicable fixes; stage your changes first
[group("lint")]
lint-fix-rs:
    {{ lint_env }} cargo clippy --locked --workspace --all-targets --fix --allow-staged

# TypeScript lints (oxlint) and the React 18 pin
[group("lint")]
lint-ts: _js
    bun run --silent lint
    jq -e '.workspaces.catalog.react | startswith("18.")' package.json >/dev/null || { printf 'the react catalog entry must stay on 18.x: Auspice runs in-process and requires React 18\n' >&2; exit 1; }

# Apply oxlint's automatic fixes
[group("lint")]
lint-fix-ts: _js
    bun run --silent lint:fix

# Custom lint libraries (dylint), denying warnings and gated against mordant-baseline.toml, then the unused public items they recorded
[group("lint")]
dylint *args:
    rm -f {{ quote(dylint_over_baseline) }}
    DYLINT_RUSTFLAGS="-A unknown_lints --deny warnings" {{ dylint_cmd }} -- {{ dylint_cargo_args }} --keep-going "$@"
    if [[ -s {{ quote(dylint_over_baseline) }} ]]; then printf 'mordant: findings over the committed baseline:\n' >&2; cat {{ quote(dylint_over_baseline) }} >&2; exit 1; fi
    cd dev/lints/dylint-custom && {{ uncached_env }} {{ pub_unused_env }} cargo run --quiet --release --locked --target-dir {{ quote(dylint_target_dir / "report") }} --bin pub-unused-report -- {{ quote(justfile_directory() / "Cargo.toml") }} {{ prepend("--exclude-crate ", public_api_crates) }}

# Apply the automatic fixes of the custom lint libraries (dylint); stage your changes first
[group("lint")]
dylint-fix:
    rm -f {{ quote(dylint_over_baseline) }}
    DYLINT_RUSTFLAGS="-A unknown_lints" {{ dylint_cmd }} --fix -- --allow-staged {{ dylint_cargo_args }}

# Accept the current mordant findings: rewrites mordant-baseline.toml, commit it afterwards
[confirm("Rewrite mordant-baseline.toml with the current findings?")]
[group("lint")]
dylint-baseline:
    DYLINT_RUSTFLAGS="-A unknown_lints" MORDANT_BASELINE_WRITE=1 {{ dylint_cmd }} -- {{ dylint_cargo_args }} --keep-going

# Unnecessary public surface (cargo-hawk), denying warnings
[group("lint")]
hawk *args:
    {{ uncached_env }} cargo +{{ hawk_toolchain }} hawk check --target-dir {{ quote(hawk_target_dir) }} {{ prepend("--exclude-crate=", public_api_crates) }} -W warnings "$@"

# Dependency bans, licenses, and sources (cargo-deny, offline)
[group("lint")]
deny:
    cargo deny --locked --all-features check bans licenses sources

# Unused Rust dependencies (cargo-shear)
[group("lint")]
shear:
    cargo shear

# TypeScript type checks of the packages, the tool configs, and the vendored lint rules
[group("lint")]
typecheck: _js
    bun run --silent typecheck

# Unused TypeScript files, exports, and dependencies (knip)
[group("lint")]
knip: _js
    bun run --silent knip
    bun run --silent knip:production

# Shell scripts in dev/ (shellcheck)
[group("lint")]
lint-shell:
    shellcheck --source-path=SCRIPTDIR $(dev/shell-files)

# Dockerfiles (hadolint)
[group("lint")]
lint-docker:
    hadolint dev/docker/*.dockerfile

# GitHub Actions workflows (actionlint)
[group("lint")]
lint-workflows:
    actionlint

# Format Rust, TypeScript, shell, TOML, and the justfile
[group("format")]
fmt: fmt-rs fmt-ts fmt-other

# Check formatting of Rust, TypeScript, shell, TOML, and the justfile
[group("format")]
fmt-check: fmt-check-rs fmt-check-ts fmt-check-other

# Format Rust (rustfmt)
[group("format")]
fmt-rs:
    cargo fmt --all

# Check formatting of Rust (rustfmt)
[group("format")]
fmt-check-rs:
    cargo fmt --all --check

# Format TypeScript, JavaScript, JSON, YAML, and CSS (oxfmt)
[group("format")]
fmt-ts: _js
    bun run --silent format

# Check formatting of TypeScript, JavaScript, JSON, YAML, and CSS (oxfmt)
[group("format")]
fmt-check-ts: _js
    bun run --silent format:check

# Format shell scripts (shfmt), TOML (taplo), and the justfile
[group("format")]
fmt-other:
    shfmt --write $(dev/shell-files)
    RUST_LOG=warn taplo fmt
    just --fmt

# Check formatting of shell scripts (shfmt), TOML (taplo), and the justfile
[group("format")]
fmt-check-other:
    shfmt --diff $(dev/shell-files)
    RUST_LOG=warn taplo fmt --check --diff
    just --fmt --check

# Regenerate the committed generated files: just gen [schemas|openapi|ts-client|napi-types|cli-docs]...
[group("generated")]
gen *groups: _js
    dev/generated generate "$@"

# Fail when a committed generated file differs from a fresh generation
[group("generated")]
generated-check: _js
    dev/generated check

# Fail when a reference fixture is missing, stale, or unclaimed in dev/registry/reference-files.toml
[group("generated")]
fixtures-check:
    dev/fixtures-check

# Start the web and API dev servers in the foreground until Ctrl-C
[group("app")]
up: _js
    dev/app up

# Probe the dev servers; exits non-zero when one is absent or built from another commit
[group("app")]
health:
    dev/app health

# Resolved ports and URLs of this checkout, and what each dev server was built from
[group("app")]
status:
    dev/app status

# Start the desktop app in development mode (in the container: TREETIME_DOCKER_X11=1)
[group("app")]
desktop: _js
    bun run --silent dev:desktop

# Production build of the web app
[group("app")]
build-web: _js
    bun run --silent build:web

# Production build of the desktop app
[group("app")]
build-desktop: _js
    bun run --silent build:desktop

# Run all benchmarks (release profile)
[group("bench")]
bench *args:
    cargo bench --locked --workspace --benches "$@"

# Sample a profile of a binary (host only, needs Docker and perf or samply): just profile treetime -- <args>
[group("bench")]
profile *args:
    dev/profile "$@"

# Explain why a crate is in the dependency tree: just why <crate> [cargo tree options]
[group("deps")]
why crate *args:
    cargo tree --locked --invert "$@"

# Update Cargo.lock within the manifest ranges to releases at least 7 days old (main checkout only)
[confirm("Update Cargo.lock?")]
[group("deps")]
deps-update *args: _main-checkout
    {{ cargo_min_age }} update "$@"
    just deps-age

# Upgrade the Rust dependency pins to the newest releases at least 7 days old (main checkout only)
[confirm("Upgrade the Rust dependency pins and rewrite Cargo.lock?")]
[group("deps")]
deps-upgrade *args: _main-checkout
    cargo upgrade --pinned --incompatible --recursive "$@"
    {{ cargo_min_age }} update
    just deps-age

# Upgrade the named JavaScript packages to their newest releases at least 7 days old (main checkout only)
[group("deps")]
deps-upgrade-ts +packages: _main-checkout
    bun update --latest "$@"

# Fail when a crate in Cargo.lock was published less than 7 days ago (network)
[group("deps")]
deps-age:
    dev/crate-age Cargo.lock dev/lints/*/Cargo.lock

# Security advisories of both dependency graphs and the crate publish age (network)
[group("deps")]
audit: audit-rs audit-ts deps-age

# Security advisories of the Rust dependencies (cargo-deny, network)
[group("deps")]
audit-rs:
    cargo deny --locked --all-features check advisories

# Security advisories of the JavaScript dependencies (bun audit, network)
[group("deps")]
audit-ts: _js
    bun audit

# Newer releases of the tools pinned in mise.toml
[group("deps")]
tools-outdated:
    mise outdated --bump

# Lock mise.lock for every supported host platform after a change to mise.toml
[group("deps")]
tools-lock:
    mise lock --platform linux-x64,linux-arm64,macos-x64,macos-arm64

# Install the JavaScript dependencies from bun.lock
_js:
    [[ -n "${TREETIME_JS_READY:-}" ]] || bun install --frozen-lockfile --silent

_main-checkout:
    test "$(git rev-parse --path-format=absolute --git-common-dir)" = "$(git rev-parse --path-format=absolute --git-dir)" || { printf 'run this recipe in the main checkout, not in a worktree\n' >&2; exit 1; }
