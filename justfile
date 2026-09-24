# TreeTime task runner: `just` lists every recipe by group.
#
# Recipes run the tools pinned in mise.toml, on the host after `mise install` or
# in the build container through `dev/docker/run just <recipe>`. Before a change
# is merged, `just check-all` must pass.

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

check_fast := "fmt-check lint typecheck lint-js"
check_full := check_fast + " dylint hawk deny shear lint-shell lint-docker lint-toml lint-workflows knip test test-js generated-check fixtures-check"

alias b := build
alias br := build-release
alias bp := build-profiling
alias r := run
alias rr := run-release
alias t := test
alias tu := test-unit
alias ti := test-integration
alias sm := smoke
alias l := lint
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

# Apply every automatic lint fix, then format; stage your changes first
[group("check")]
fix: lint-fix fmt

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

# Rust tests (nextest); arguments are nextest filters and options
[group("test")]
test *args:
    cargo nextest run --locked --workspace "$@"

# Rust unit tests only
[group("test")]
test-unit *args:
    cargo nextest run --locked --workspace --lib "$@"

# Rust integration tests only
[group("test")]
test-integration *args:
    cargo nextest run --locked --workspace --test '*' "$@"

# List Rust tests without running them
[group("test")]
test-list *args:
    cargo nextest list --locked --workspace "$@"

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

# TypeScript tests (vitest) and the custom oxlint rule tests
[group("test")]
test-js: _js
    bun run --silent test
    node --test "dev/lints/oxlint/__tests__/test_*.ts" "dev/lints/oxlint-anti-slop/**/*.test.ts"

# Coverage report, never a gate: just coverage [rs|ts|all]
[group("report")]
coverage target="all": _js
    case {{ quote(target) }} in rs | ts | all) ;; *) printf 'unknown coverage target: %s (rs, ts, all)\n' {{ quote(target) }} >&2; exit 1 ;; esac
    if [[ {{ quote(target) }} != ts ]]; then {{ uncached_env }} cargo llvm-cov nextest --locked --workspace --html --output-dir tmp/coverage/rs; printf 'Rust coverage: tmp/coverage/rs/html/index.html\n'; fi
    if [[ {{ quote(target) }} != rs ]]; then bun run --silent coverage; printf 'TypeScript coverage: tmp/coverage/ts/index.html\n'; fi

# Mutation testing of the code changed since the fork point (or since `rust`), never a gate
[group("report")]
mutants *args:
    dev/mutants "$@"

# Copy-paste duplication across Rust and TypeScript (jscpd), never a gate
[group("report")]
duplication *args:
    jscpd --config .jscpd.json "$@"
# Explain why a crate is in the dependency tree: just why <crate> [cargo tree options]
why crate *args:
    cargo tree --locked --invert "$@"

# Inventory of lint suppressions and other review-sensitive settings
[group("report")]
review-suppressions:
    dev/review-suppressions

# Clippy over all targets, denying warnings
[group("lint")]
lint *args:
    {{ lint_env }} cargo clippy --locked --workspace --all-targets "$@" -- --deny warnings

# Apply clippy's machine-applicable fixes, then the dylint and oxlint fixes; stage your changes first
[group("lint")]
lint-fix: _js
    {{ lint_env }} cargo clippy --locked --workspace --all-targets --fix --allow-staged
    rm -f {{ quote(dylint_target_dir / "mordant/over-baseline.txt") }}
    {{ uncached_env }} CARGO_TARGET_DIR={{ quote(dylint_target_dir) }} DYLINT_RUSTFLAGS="-A unknown_lints" TREETIME_LINTS_PUB_UNUSED_DIR={{ quote(dylint_target_dir / "pub-unused") }} cargo dylint --quiet --all --fix -- --allow-staged --quiet --locked --workspace --all-targets
    bun run --silent lint:fix

# Custom lint libraries (dylint), denying warnings and gated against mordant-baseline.toml
[group("lint")]
dylint *args:
    rm -f {{ quote(dylint_target_dir / "mordant/over-baseline.txt") }}
    {{ uncached_env }} CARGO_TARGET_DIR={{ quote(dylint_target_dir) }} DYLINT_RUSTFLAGS="-A unknown_lints --deny warnings" TREETIME_LINTS_PUB_UNUSED_DIR={{ quote(dylint_target_dir / "pub-unused") }} cargo dylint --quiet --all -- --quiet --locked --workspace --all-targets --keep-going "$@"
    if [[ -s {{ quote(dylint_target_dir / "mordant/over-baseline.txt") }} ]]; then printf 'mordant: findings over the committed baseline:\n' >&2; cat {{ quote(dylint_target_dir / "mordant/over-baseline.txt") }} >&2; exit 1; fi
    cd dev/lints/dylint-custom && {{ uncached_env }} TREETIME_LINTS_PUB_UNUSED_DIR={{ quote(dylint_target_dir / "pub-unused") }} cargo run --quiet --release --locked --target-dir {{ quote(dylint_target_dir / "report") }} --bin pub-unused-report -- {{ quote(justfile_directory() / "Cargo.toml") }} {{ prepend("--exclude-crate ", public_api_crates) }}
    cd dev/lints/dylint-custom && {{ uncached_env }} cargo test --quiet --release --locked --lib --bins --target-dir {{ quote(dylint_target_dir / "report") }}
    cd dev/lints/dylint-trailofbits && {{ uncached_env }} cargo test --quiet --locked --workspace --lib --target-dir {{ quote(dylint_target_dir / "trailofbits-test") }}

# Accept the current mordant findings: rewrites mordant-baseline.toml, commit it afterwards
[confirm("Rewrite mordant-baseline.toml with the current findings?")]
[group("lint")]
dylint-baseline:
    {{ uncached_env }} CARGO_TARGET_DIR={{ quote(dylint_target_dir) }} DYLINT_RUSTFLAGS="-A unknown_lints" MORDANT_BASELINE_WRITE=1 cargo dylint --quiet --all -- --quiet --keep-going --locked --workspace --all-targets

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

# TypeScript lints (oxlint) and the React 18 pin
[group("lint")]
lint-js: _js
    bun run --silent lint
    jq -e '.workspaces.catalog.react | startswith("18.")' package.json >/dev/null || { printf 'the react catalog entry must stay on 18.x: Auspice runs in-process and requires React 18\n' >&2; exit 1; }

# TypeScript type checks of the packages, the tool configs, and the vendored lint rules
[group("lint")]
typecheck: _js
    bun run --silent typecheck

# Unused TypeScript files, exports, and dependencies (knip)
[group("lint")]
knip: _js
    bun run --silent knip
    bun run --silent knip:production

# Shell scripts in dev/ (shellcheck, shfmt)
[group("lint")]
lint-shell:
    shellcheck --source-path=SCRIPTDIR $(dev/shell-files)
    shfmt --diff $(dev/shell-files)

# Dockerfiles (hadolint)
[group("lint")]
lint-docker:
    hadolint dev/docker/*.dockerfile

# TOML formatting (taplo)
[group("lint")]
lint-toml:
    RUST_LOG=warn taplo fmt --check --diff

# GitHub Actions workflows (actionlint)
[group("lint")]
lint-workflows:
    actionlint

# Format Rust, TypeScript, shell, TOML, and the justfile
[group("check")]
fmt: _js
    cargo fmt --all
    bun run --silent format
    shfmt --write $(dev/shell-files)
    RUST_LOG=warn taplo fmt
    just --fmt

# Check formatting of Rust, TypeScript, shell, TOML, and the justfile
[group("check")]
fmt-check: _js
    cargo fmt --all --check
    bun run --silent format:check
    shfmt --diff $(dev/shell-files)
    RUST_LOG=warn taplo fmt --check
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
deps-upgrade-js +packages: _main-checkout
    bun update --latest "$@"

# Fail when a crate in Cargo.lock was published less than 7 days ago (network)
[group("deps")]
deps-age:
    dev/crate-age Cargo.lock dev/lints/*/Cargo.lock

# Security advisories of both dependency graphs and the crate publish age (network)
[group("deps")]
audit: _js
    cargo deny --locked --all-features check advisories
    bun audit
    just deps-age

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
