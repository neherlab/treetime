# TreeTime task runner. Run through the container boundary, e.g.
#   ./dev/docker/run just check
# `just --list` shows every task grouped with a one-line description.
# Tool versions come from mise (mise.toml / mise.lock).

set shell := ["bash", "-euo", "pipefail", "-c"]
set positional-arguments := true
export CARGO_TERM_QUIET := "true"
export KACHE_PROGRESS := "off"

project_dir := justfile_directory()

# Per-kind cargo target directories. The container bind-mounts .build, so these
# match the paths dev/docker/run mounts and keep the build cache warm.
build_dir := project_dir / ".build/docker"
test_dir := build_dir / "test"

# Isolated target dir for the dylint driver: the pinned nightly toolchain must
# not overwrite the stable build cache.
dylint_dir := build_dir / "dylint"

# Records the pub_unused_in_workspace lint writes per crate for pub-unused-report.
pub_unused_dir := dylint_dir / "treetime_lints/pub_unused"

# Show the grouped task list (default).
default:
    @just --list --list-heading $'TreeTime tasks (run via ./dev/docker/run just <task>):\n'

# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

alias b := build
alias br := build-release
alias bp := build-profiling

# Build the workspace (dev profile)
[group('build')]
build *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{build_dir}}' RUSTFLAGS="$(rustflags_build)"
    kache_use build
    nicely cargo -q build --locked "$@"

# Build the workspace (release profile); copies the named binary to .out
[group('build')]
build-release *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{build_dir}}' RUSTFLAGS="$(rustflags_build)"
    kache_use build
    bin="${1:-}"; [[ $# -gt 0 ]] && shift || true
    nicely cargo -q build --release --locked ${bin:+--bin="${bin}"} "$@"
    [[ -n "${bin}" ]] && copy_bin_to_out "${bin}" || true

# Build the workspace (profiling profile: release + debug symbols)
[group('build')]
build-profiling *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{build_dir}}' RUSTFLAGS="$(rustflags_build)"
    kache_use build
    bin="${1:-}"; [[ $# -gt 0 ]] && shift || true
    nicely cargo -q build --profile=profiling --locked ${bin:+--bin="${bin}"} "$@"

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

alias r := run
alias rr := run-release
alias E := example
alias Er := example-release

# Run a binary (dev profile): just run <bin> [args...]
[group('run')]
run bin *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{build_dir}}' RUSTFLAGS="$(rustflags_build)"
    kache_use build
    nicely cargo -q run --locked --bin "$@"

# Run a binary (release profile): just run-release <bin> [args...]
[group('run')]
run-release bin *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{build_dir}}' RUSTFLAGS="$(rustflags_build)"
    kache_use build
    nicely cargo -q run --release --locked --bin "$@"

# Run an example (dev profile): just example <name> [args...]
[group('run')]
example name *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{build_dir}}' RUSTFLAGS="$(rustflags_build)"
    kache_use build
    nicely cargo -q run --locked --example "$@"

# Run an example (release profile): just example-release <name> [args...]
[group('run')]
example-release name *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{build_dir}}' RUSTFLAGS="$(rustflags_build)"
    kache_use build
    nicely cargo -q run --release --locked --example "$@"

# ---------------------------------------------------------------------------
# Test
# ---------------------------------------------------------------------------

alias t := test-all
alias tu := test-unit
alias ti := test-integration
alias tl := test-list

# Run all tests (nextest); fails when a filter matches nothing
[group('test')]
test-all *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{test_dir}}' RUSTFLAGS="$(rustflags_test)"
    kache_use test
    nextest_guard "$@"
    nicely cargo -q nextest run --locked --success-output=immediate --workspace --cargo-quiet --no-fail-fast --hide-progress-bar --color=always "$@" 2>&1 \
      | '{{project_dir}}/dev/prettytest'

# Run unit (lib) tests only; fails when a filter matches nothing
[group('test')]
test-unit *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{test_dir}}' RUSTFLAGS="$(rustflags_test)"
    kache_use test
    nextest_guard --lib "$@"
    nicely cargo -q nextest run --lib --locked --success-output=immediate --workspace --cargo-quiet --no-fail-fast --hide-progress-bar --color=always "$@" 2>&1 \
      | '{{project_dir}}/dev/prettytest'

# Run integration tests only; fails when a filter matches nothing
[group('test')]
test-integration *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{test_dir}}' RUSTFLAGS="$(rustflags_test)"
    kache_use test
    nextest_guard --test='*' "$@"
    nicely cargo -q nextest run --test='*' --locked --success-output=immediate --workspace --cargo-quiet --no-fail-fast --hide-progress-bar --color=always "$@" 2>&1 \
      | '{{project_dir}}/dev/prettytest'

# List tests without running them
[group('test')]
test-list *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{test_dir}}' RUSTFLAGS="$(rustflags_test)"
    kache_use test
    pretty_args=(); cargo_args=()
    while (( $# > 0 )); do
      case "${1}" in
        --summary|--summary=*|--filter=*) pretty_args+=("${1}"); shift ;;
        --) shift; break ;;
        -*) cargo_args+=("${1}"); shift ;;
        *) pretty_args+=("${1}"); shift ;;
      esac
    done
    cargo_args+=("$@")
    cargo -q nextest list --locked --workspace --cargo-quiet "${cargo_args[@]}" \
      | '{{project_dir}}/dev/prettytestlist' "${pretty_args[@]}"

# ---------------------------------------------------------------------------
# Coverage
# ---------------------------------------------------------------------------

alias cov := coverage
alias cl := coverage-lcov

# HTML coverage report under tmp/coverage
[group('test')]
coverage *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{test_dir}}' RUSTFLAGS="$(rustflags_test)"
    kache_use test
    cov_dir='{{project_dir}}/tmp/coverage'; mkdir -p "${cov_dir}"
    nicely cargo llvm-cov nextest --html --output-dir "${cov_dir}" --ignore-run-fail --locked --workspace --hide-progress-bar --color=always "$@"
    printf 'Coverage report: %s/index.html\n' "${cov_dir}"

# LCOV coverage report under tmp/coverage
[group('test')]
coverage-lcov *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{test_dir}}' RUSTFLAGS="$(rustflags_test)"
    kache_use test
    cov_dir='{{project_dir}}/tmp/coverage'; mkdir -p "${cov_dir}"
    nicely cargo llvm-cov nextest --lcov --output-path "${cov_dir}/lcov.info" --ignore-run-fail --locked --workspace --hide-progress-bar --color=always "$@"
    printf 'LCOV: %s/lcov.info\n' "${cov_dir}"

alias cu := coverage-uncovered
alias cj := coverage-js

# List Rust source files that have uncovered lines (reports, never gates)
[group('test')]
coverage-uncovered *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{test_dir}}' RUSTFLAGS="$(rustflags_test)"
    kache_use test
    cov_dir='{{project_dir}}/tmp/coverage'; mkdir -p "${cov_dir}"
    summary="${cov_dir}/summary.json"
    nicely cargo llvm-cov nextest --json --summary-only --output-path "${summary}" --ignore-run-fail --locked --workspace --hide-progress-bar "$@"
    printf '\nFiles with uncovered lines (uncovered/total lines):\n'
    jq -r '.data[0].files[] | select(.summary.lines.count > .summary.lines.covered) | "\(.summary.lines.count - .summary.lines.covered)\t\(.summary.lines.count)\t\(.filename)"' "${summary}" \
      | sort -rn | awk -F'\t' '{printf "  %6s/%-6s %s\n", $1, $2, $3}'
    uncovered_files="$(jq -r '[.data[0].files[] | select(.summary.lines.count > .summary.lines.covered)] | length' "${summary}")"
    total_files="$(jq -r '.data[0].files | length' "${summary}")"
    printf '\n%s of %s files have uncovered lines. Full summary: %s\n' "${uncovered_files}" "${total_files}" "${summary}"

# Report JavaScript/TypeScript coverage (vitest v8; reports, never gates)
[group('test')]
coverage-js *args:
    #!/usr/bin/env bash
    set -euo pipefail
    cd '{{project_dir}}'
    if [[ ! -x node_modules/.bin/vitest ]]; then
      printf 'vitest is not installed; run `just setup` and set up the JS test runner first. No JS coverage to report.\n' >&2
      exit 0
    fi
    bun run coverage "$@"

# ---------------------------------------------------------------------------
# Mutation testing
# ---------------------------------------------------------------------------

alias mut := mutants
alias mutf := mutants-full

# Mutation testing on the fork-point diff (cargo-mutants via nextest)
[group('test')]
mutants *args:
    @just _mutants diff "$@"

# Mutation testing across the whole workspace (cargo-mutants via nextest)
[group('test')]
mutants-full *args:
    @just _mutants full "$@"

# Shared mutation driver: validate exclusions, resolve scope, run cargo-mutants.
_mutants scope *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    cd '{{project_dir}}'
    # positional-arguments passes scope as $1 too; drop it so "$@" is user args.
    shift
    mkdir -p tmp/mutants
    native='.cargo/mutants.toml'
    excl_file='.cargo/mutants-exclusions.toml'

    # Guard: raw exclude/examine keys in the native config would bypass the
    # reason requirement, so reject them.
    if yq -p toml -o json '.' "${native}" \
      | jq -e 'has("exclude_re") or has("exclude_globs") or has("examine_re") or has("examine_globs")' >/dev/null; then
      printf 'error: %s must not declare exclude_*/examine_* keys; declare exclusions with a reason in %s\n' "${native}" "${excl_file}" >&2
      exit 1
    fi

    # Validate every exclusion carries a reason and exactly one of file/name,
    # and turn each into a cargo-mutants flag.
    excl_args=()
    count="$(yq -p toml -oy '.exclude | length' "${excl_file}")"
    for ((i = 0; i < count; i++)); do
      reason="$(yq -p toml -oy ".exclude[${i}].reason // \"\"" "${excl_file}")"
      file="$(yq -p toml -oy ".exclude[${i}].file // \"\"" "${excl_file}")"
      name="$(yq -p toml -oy ".exclude[${i}].name // \"\"" "${excl_file}")"
      if [[ -z "${reason}" || "${reason}" == "null" ]]; then
        printf 'error: mutation exclusion #%s in %s has no reason; a reason is required\n' "${i}" "${excl_file}" >&2
        exit 1
      fi
      if [[ -n "${file}" && -n "${name}" ]] || [[ -z "${file}" && -z "${name}" ]]; then
        printf 'error: mutation exclusion #%s in %s must set exactly one of file or name\n' "${i}" "${excl_file}" >&2
        exit 1
      fi
      [[ -n "${file}" ]] && excl_args+=(--exclude "${file}")
      [[ -n "${name}" ]] && excl_args+=(--exclude-re "${name}")
    done

    scope_args=()
    if [[ '{{scope}}' == "diff" ]]; then
      branch="$(git branch --show-current)"
      # The fork point is what create-worktree recorded; never the base branch.
      if command -v create-worktree >/dev/null 2>&1; then
        fork="$(create-worktree --fork-point 2>/dev/null || true)"
      fi
      fork="${fork:-$(git config "branch.${branch}.fork-point" 2>/dev/null || true)}"
      if [[ -z "${fork}" ]]; then
        printf 'error: no fork point for %s (branch.%s.fork-point unset and create-worktree unavailable)\n' "${branch}" "${branch}" >&2
        exit 1
      fi
      diff_file="$(mktemp)"
      trap 'rm -f "${diff_file}"' EXIT
      git diff "${fork}..HEAD" >"${diff_file}"
      printf 'Mutating the diff against fork point %s\n' "${fork}" >&2
      scope_args=(--in-diff "${diff_file}")
    fi

    # cargo-mutants manages its own per-mutant build directories, so do not pin
    # CARGO_TARGET_DIR (that would make parallel mutant builds overwrite each
    # other). RUSTFLAGS still carries the openblas link flags the tests need.
    export RUSTFLAGS="$(rustflags_test)"
    nicely cargo mutants --test-tool nextest --colors always "${scope_args[@]}" "${excl_args[@]}" "$@"

# ---------------------------------------------------------------------------
# Duplication
# ---------------------------------------------------------------------------

alias dup := duplication

# Report copy-paste duplication across Rust and TypeScript (jscpd, reports only)
[group('test')]
duplication *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    cd '{{project_dir}}'
    mkdir -p tmp/jscpd
    nicely jscpd --config .jscpd.json "$@"

# ---------------------------------------------------------------------------
# Lint and format (Rust)
# ---------------------------------------------------------------------------

alias l := lint
alias lf := lint-fix
alias lc := lint-ci
alias lD := lint-deps
alias dy := dylint-all
alias dyf := dylint-all-fix
alias f := format
alias fc := format-check
alias q := quality

# Clippy over all targets
[group('lint')]
clippy *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{build_dir}}' RUSTFLAGS="$(rustflags_build)"
    kache_use clippy
    nicely cargo -q clippy -q --all-targets --all --locked "$@"

# Clippy with autofix (stage changes first)
[group('lint')]
clippy-fix *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{build_dir}}' RUSTFLAGS="$(rustflags_build)"
    kache_use clippy
    vcs_flag="--allow-staged"; [[ -f '{{project_dir}}/.git' ]] && vcs_flag="--allow-no-vcs"
    nicely cargo -q clippy -q --all-targets --all --fix "${vcs_flag}" --locked "$@"

# Lint: clippy
[group('lint')]
lint *args:
    #!/usr/bin/env bash
    set -euo pipefail
    just clippy "$@"

# Lint with autofix
[group('lint')]
lint-fix *args:
    #!/usr/bin/env bash
    set -euo pipefail
    just clippy-fix "$@"

# Lint: clippy + all dylint libraries
[group('lint')]
lint-extra *args:
    #!/usr/bin/env bash
    set -euo pipefail
    just clippy "$@"
    just dylint-all

# Lint with autofix: clippy + every machine-applicable Dylint fix
[group('lint')]
lint-extra-fix *args:
    #!/usr/bin/env bash
    set -euo pipefail
    just clippy-fix "$@"
    just dylint-all-fix

# Lint: every Rust checker (clippy, dylint, hawk, deny, shear)
[group('lint')]
lint-all *args:
    #!/usr/bin/env bash
    set -euo pipefail
    just clippy "$@"
    just dylint-all
    just hawk -W warnings
    just lint-deps
    if [[ -f '{{project_dir}}/deny.toml' ]] && command -v cargo-deny >/dev/null 2>&1; then
      cargo deny --locked check bans licenses sources
    fi

# Lint with autofix: every Rust checker that supports it
[group('lint')]
lint-all-fix *args:
    #!/usr/bin/env bash
    set -euo pipefail
    just clippy-fix "$@"
    just dylint-all-fix

# Clippy denying warnings (CI parity)
[group('lint')]
lint-ci *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{build_dir}}' RUSTFLAGS="$(rustflags_build)"
    kache_use lint-ci
    nicely cargo -q clippy -q --all-targets --all --locked -- -Dwarnings "$@"

# Report unused workspace dependencies (cargo-shear)
[group('lint')]
lint-deps *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    nicely cargo shear "$@"

# Run every Dylint library declared in workspace metadata
[group('lint')]
dylint-all *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    # No compiler cache: a cache hit skips the lint passes, so mordant's
    # baseline check and the pub_unused_in_workspace records would not be written.
    unset RUSTFLAGS RUSTC_WRAPPER
    export CARGO_TARGET_DIR='{{dylint_dir}}' DYLINT_RUSTFLAGS="-A unknown_lints" RUST_BACKTRACE=0 CARGO_INCREMENTAL=0
    export TREETIME_LINTS_PUB_UNUSED_DIR='{{pub_unused_dir}}'
    rm -f '{{dylint_dir}}/mordant/over-baseline.txt'
    nicely cargo dylint --quiet --all -- --quiet --locked --workspace --all-targets --keep-going "$@"
    status=0
    just _pub-unused-report || status=1
    if [[ -s '{{dylint_dir}}/mordant/over-baseline.txt' ]]; then
      printf 'mordant: findings over the committed baseline:\n' >&2
      cat '{{dylint_dir}}/mordant/over-baseline.txt' >&2
      status=1
    fi
    exit "${status}"

# Run one lint from the vendored Trail of Bits library
[group('lint')]
dylint-trailofbits lint *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    unset RUSTFLAGS
    lint="${1}"
    shift
    allow_others=()
    for manifest in '{{project_dir}}'/dev/lints/dylint-trailofbits/*/Cargo.toml; do
      other="$(basename "$(dirname "${manifest}")")"
      if [[ "${other}" != "${lint}" ]]; then
        allow_others+=(-A "${other}")
      fi
    done
    export CARGO_TARGET_DIR='{{dylint_dir}}' DYLINT_RUSTFLAGS="-A unknown_lints ${allow_others[*]} -W ${lint}" RUST_BACKTRACE=0 CARGO_INCREMENTAL=0
    kache_use dylint
    nicely cargo dylint --quiet --lib trailofbits -- --quiet --locked --workspace --all-targets "$@"

# Apply every machine-applicable Dylint fix
[group('lint')]
dylint-all-fix *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    # No compiler cache, as in dylint-all.
    unset RUSTFLAGS RUSTC_WRAPPER
    export CARGO_TARGET_DIR='{{dylint_dir}}' DYLINT_RUSTFLAGS="-A unknown_lints" RUST_BACKTRACE=0 CARGO_INCREMENTAL=0
    export TREETIME_LINTS_PUB_UNUSED_DIR='{{pub_unused_dir}}'
    rm -f '{{dylint_dir}}/mordant/over-baseline.txt'
    vcs_flag="--allow-staged"; [[ -f '{{project_dir}}/.git' ]] && vcs_flag="--allow-no-vcs"
    nicely cargo dylint --quiet --all --fix -- "${vcs_flag}" --quiet --locked --workspace --all-targets --keep-going "$@"
    status=0
    just _pub-unused-report || status=1
    if [[ -s '{{dylint_dir}}/mordant/over-baseline.txt' ]]; then
      printf 'mordant: findings over the committed baseline:\n' >&2
      cat '{{dylint_dir}}/mordant/over-baseline.txt' >&2
      status=1
    fi
    exit "${status}"

# Run the tests of the custom Dylint library's report binaries
[group('lint')]
dylint-custom-test *args:
    #!/usr/bin/env bash
    set -euo pipefail
    unset RUSTFLAGS RUSTC_WRAPPER CARGO_TARGET_DIR
    pushd '{{project_dir}}/dev/lints/dylint-custom' >/dev/null
    cargo test --quiet --release --locked --target-dir '{{dylint_dir}}/pub-unused-report' --bins "$@"
    popd >/dev/null

# Run the UI tests of the vendored Trail of Bits Dylint library
[group('lint')]
dylint-trailofbits-test *args:
    #!/usr/bin/env bash
    set -euo pipefail
    unset RUSTFLAGS RUSTC_WRAPPER
    export CARGO_TARGET_DIR='{{dylint_dir}}/trailofbits-test'
    pushd '{{project_dir}}/dev/lints/dylint-trailofbits' >/dev/null
    cargo test --quiet --locked --workspace --lib "$@"
    popd >/dev/null

# Report public items no workspace crate uses, from the records the
# pub_unused_in_workspace lint wrote during the last dylint run
_pub-unused-report:
    #!/usr/bin/env bash
    set -euo pipefail
    unset RUSTFLAGS RUSTC_WRAPPER CARGO_TARGET_DIR
    export TREETIME_LINTS_PUB_UNUSED_DIR='{{pub_unused_dir}}'
    pushd '{{project_dir}}/dev/lints/dylint-custom' >/dev/null
    cargo run --quiet --release --locked --target-dir '{{dylint_dir}}/pub-unused-report' --bin pub-unused-report -- '{{project_dir}}/Cargo.toml'
    popd >/dev/null

# Report unnecessary public surface across the workspace (cargo-hawk)
[group('lint')]
hawk *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    unset RUSTFLAGS CARGO_TARGET_DIR
    # cargo-hawk ships a compiler driver tied to a specific rustc; run it through
    # that toolchain (installed by dev/docker/files/install-hawk, version shared
    # via dev/docker/files/hawk-toolchain).
    toolchain="$(cat '{{project_dir}}/dev/docker/files/hawk-toolchain')"
    # Library crates whose public API is an external boundary: the shared utility
    # and file-format libraries publish a complete API for reuse, including
    # operations no current caller needs, and the Node addon's `#[napi]` surface
    # is consumed by JavaScript, not by a Rust target. hawk.toml selects only
    # modules and files, so whole crates are excluded here.
    excluded_crates=(
      treetime_utils
      util_newick
      util_phyloxml
      util_augur_node_data_json
      util_usher_mat
      app_napi
    )
    exclude_flags=()
    for crate in "${excluded_crates[@]}"; do
      exclude_flags+=(--exclude-crate="${crate}")
    done
    nicely cargo "+${toolchain}" hawk check --target-dir '{{build_dir}}/hawk' "${exclude_flags[@]}" "$@"

# Regenerate the committed mordant baseline (mordant-baseline.toml)
[group('lint')]
dylint-mordant-baseline *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    # Match the mordant gate: fixed (empty) RUSTFLAGS so the seeded baseline and
    # the check run analyze the workspace identically.
    unset RUSTFLAGS RUSTC_WRAPPER
    export CARGO_TARGET_DIR='{{dylint_dir}}' DYLINT_RUSTFLAGS="-A unknown_lints" MORDANT_BASELINE_WRITE=1 CARGO_INCREMENTAL=0
    nicely cargo dylint --quiet --all -- --quiet --keep-going --locked --workspace --all-targets "$@"
    printf 'Regenerated mordant-baseline.toml\n'

# Lint levels, allow/expect lists, mutation exclusions, ignored tests, and float
# tolerances are changed rarely and reviewed separately from ordinary code.

# Inventory lint suppressions and other review-sensitive surfaces (read-only)
[group('lint')]
review-suppressions *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    cd '{{project_dir}}'
    report() { printf '\n=== %s ===\n' "$1"; shift; "$@" || true; }
    # Explicit source roots keep ripgrep out of the build and vendor trees.
    report "rust #[allow]/#[expect] attributes" \
      rg -n -g '*.rs' -g '!**/generated/**' '#!?\[(allow|expect)\(' packages
    report "clippy/rustc allow entries in manifests" \
      rg -n -g 'Cargo.toml' -e '= "allow"' -e 'level *= *"allow"' Cargo.toml packages
    report "ignored or skipped rust tests" \
      rg -n -g '*.rs' '#\[ignore' packages
    report "cargo-mutants exclusions" \
      rg -n -g 'mutants.toml' -e 'exclude' -e 'skip' .cargo .config
    report "float comparison tolerances" \
      rg -n -g '*.rs' 'epsilon *= *1e-|max_ulps *= *' packages
    report "typescript lint suppressions" \
      rg -n -g '*.ts' -g '*.tsx' -g '!**/generated/**' 'oxlint-disable|eslint-disable|@ts-(expect-error|ignore)' packages
    printf '\nReview these separately from ordinary code changes.\n'

# Format Rust code
[group('lint')]
format *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    nicely cargo -q fmt "$@"

# Check Rust formatting (read-only)
[group('lint')]
format-check *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    nicely cargo -q fmt --check "$@"

# lint-fix, then format, then test
[group('lint')]
quality *args:
    #!/usr/bin/env bash
    set -euo pipefail
    just lint-fix "$@"
    just format
    just test-all "$@"

# ---------------------------------------------------------------------------
# Read-only gates
# ---------------------------------------------------------------------------

# Fast read-only checks: rust format, clippy, TypeScript, oxlint
[group('gate')]
check: (_check "fast")

# Full read-only checks: fast set plus dylint, deny, unused code, freshness, tests
[group('gate')]
check-all: (_check "full")

_check mode:
    #!/usr/bin/env bash
    set -uo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    cd '{{project_dir}}'
    names=(); results=()
    record() { names+=("$1"); results+=("$2"); printf '  -> %s: %s\n' "$1" "$2" >&2; }
    run_check() {
      local name="$1"; shift
      printf '\n=== %s ===\n' "$name" >&2
      if "$@"; then record "$name" PASS; else record "$name" "FAIL (exit $?)"; fi
    }
    skip() { record "$1" "SKIP ($2)"; }

    have_bun_project() { command -v bun >/dev/null 2>&1 && [[ -x '{{project_dir}}/node_modules/.bin/turbo' ]]; }

    kache_base="${KACHE_CACHE_DIR:-}"
    export CARGO_TARGET_DIR='{{build_dir}}' RUSTFLAGS="$(rustflags_build)"
    export KACHE_CACHE_DIR="${kache_base}"; kache_use clippy

    run_check "rust-format" cargo -q fmt --all --check
    if [[ "{{mode}}" == "full" ]]; then
      run_check "rust-lints" just lint-all
    else
      run_check "rust-clippy" just clippy
    fi
    if have_bun_project; then
      run_check "typescript" bash -c "cd '{{project_dir}}' && bun run typecheck"
      run_check "typescript-config" bash -c "cd '{{project_dir}}' && bun run typecheck:tools"
      run_check "typescript-vendor" bash -c "cd '{{project_dir}}' && bun run typecheck:vendor"
      run_check "oxlint" bash -c "cd '{{project_dir}}' && bun run lint"
    else
      skip "typescript" "no node_modules; run just setup"
      skip "typescript-config" "no node_modules; run just setup"
      skip "typescript-vendor" "no node_modules; run just setup"
      skip "oxlint" "no node_modules; run just setup"
    fi

    if [[ "{{mode}}" == "full" ]]; then
      if command -v shellcheck >/dev/null 2>&1; then
        mapfile -t sh_files < <(dev_shell_files '{{project_dir}}')
        run_check "shellcheck" shellcheck --source-path='{{project_dir}}/dev/cross:{{project_dir}}/dev/docker:{{project_dir}}/dev/lib:{{project_dir}}/dev' "${sh_files[@]}"
      else
        skip "shellcheck" "shellcheck not installed"
      fi
      if command -v shfmt >/dev/null 2>&1; then
        mapfile -t sh_files < <(dev_shell_files '{{project_dir}}')
        run_check "shfmt" shfmt --diff "${sh_files[@]}"
      else
        skip "shfmt" "shfmt not installed"
      fi
      if command -v taplo >/dev/null 2>&1; then
        run_check "toml-format" taplo fmt --check --diff
      else
        skip "toml-format" "taplo not installed"
      fi
      if command -v hadolint >/dev/null 2>&1; then
        shopt -s nullglob
        run_check "hadolint" hadolint '{{project_dir}}'/dev/docker/*.dockerfile
      else
        skip "hadolint" "hadolint not installed"
      fi
      if have_bun_project && grep -q '"knip"' package.json 2>/dev/null; then
        run_check "knip" bash -c "cd '{{project_dir}}' && bun run knip"
        run_check "knip-production" bash -c "cd '{{project_dir}}' && bun run knip:production"
      else
        skip "knip" "not configured"
        skip "knip-production" "not configured"
      fi
      run_check "generated-freshness" bash -c "just generated-check"
      export CARGO_TARGET_DIR='{{test_dir}}' RUSTFLAGS="$(rustflags_test)"
      export KACHE_CACHE_DIR="${kache_base}"; kache_use test
      run_check "tests" cargo -q nextest run --locked --workspace --cargo-quiet --no-fail-fast --hide-progress-bar
      if have_bun_project; then
        run_check "js-tests" bash -c "cd '{{project_dir}}' && bun run test"
        run_check "oxlint-rules" bash -c "cd '{{project_dir}}' && node --test \"dev/lints/oxlint/__tests__/test_*.ts\" \"dev/lints/oxlint-anti-slop/**/*.test.ts\""
      else
        skip "js-tests" "no node_modules"
        skip "oxlint-rules" "no node_modules"
      fi
      run_check "dylint-custom-tests" just dylint-custom-test
      run_check "dylint-trailofbits-tests" just dylint-trailofbits-test
    fi

    printf '\n===== %s check summary =====\n' "{{mode}}"
    fail=0
    for i in "${!names[@]}"; do
      printf '  %-22s %s\n' "${names[$i]}" "${results[$i]}"
      [[ "${results[$i]}" == FAIL* ]] && fail=$((fail + 1))
    done
    printf '\n'
    if (( fail > 0 )); then
      printf '%d check(s) reported FAIL (see summary above). Skipped checks do not count as failures.\n' "${fail}" >&2
      exit 1
    fi
    printf 'All checks that ran passed.\n'

# ---------------------------------------------------------------------------
# Shell scripts (dev/)
# ---------------------------------------------------------------------------

# Lint dev/ shell scripts (shellcheck)
[group('shell')]
sh-lint *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    mapfile -t files < <(dev_shell_files '{{project_dir}}')
    nicely shellcheck --source-path='{{project_dir}}/dev/cross:{{project_dir}}/dev/docker:{{project_dir}}/dev/lib:{{project_dir}}/dev' "$@" "${files[@]}"

# Check dev/ shell script formatting (shfmt, read-only)
[group('shell')]
sh-format-check *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    mapfile -t files < <(dev_shell_files '{{project_dir}}')
    nicely shfmt --diff "$@" "${files[@]}"

# Format dev/ shell scripts in place (shfmt)
[group('shell')]
sh-format *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    mapfile -t files < <(dev_shell_files '{{project_dir}}')
    nicely shfmt --write "$@" "${files[@]}"

# ---------------------------------------------------------------------------
# TOML and Dockerfiles
# ---------------------------------------------------------------------------

# Format all TOML and sort dependency tables (taplo)
[group('config')]
toml-format *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    nicely taplo fmt "$@"

# Check TOML formatting and dependency ordering (read-only)
[group('config')]
toml-check *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    nicely taplo fmt --check --diff "$@"

# Lint Dockerfiles (hadolint)
[group('config')]
docker-lint *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    shopt -s nullglob
    files=('{{project_dir}}'/dev/docker/*.dockerfile)
    nicely hadolint "$@" "${files[@]}"

# ---------------------------------------------------------------------------
# Benchmarks and profiling
# ---------------------------------------------------------------------------

alias B := bench-all
alias Bd := bench-debug

# Run all benchmarks (release)
[group('bench')]
bench-all *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{build_dir}}' RUSTFLAGS="$(rustflags_test)"
    kache_use build
    nicely cargo -q bench --workspace --benches --locked "$@"

# Run all benchmarks (profiling profile)
[group('bench')]
bench-debug *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{build_dir}}' RUSTFLAGS="$(rustflags_test)"
    kache_use build
    nicely cargo -q bench --workspace --benches --profile=profiling --locked "$@"

# ---------------------------------------------------------------------------
# Docs and inspection
# ---------------------------------------------------------------------------

alias U := upgrade-deps
alias w := why
alias D := docs
alias L := list

# Upgrade Rust dependencies (main checkout only)
[group('deps')]
upgrade-deps *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    require_main_checkout '{{project_dir}}' "dependency upgrades"
    nicely cargo -q upgrade --pinned --incompatible --verbose --recursive "$@"

# Update Cargo.lock to the newest versions the manifests allow (main checkout only): just update-lock [-p <crate>]
[group('deps')]
update-lock *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    require_main_checkout '{{project_dir}}' "lockfile updates"
    nicely cargo update "$@"

# Explain why a crate is in the dependency tree: just why <crate>
[group('deps')]
why *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    nicely cargo -q tree --locked -i "$@"

# Security audit over the network: cargo-deny advisories (fetches the RustSec database) and bun audit
[group('deps')]
audit *args:
    #!/usr/bin/env bash
    set -uo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    cd '{{project_dir}}'
    rc=0
    printf '\n=== cargo-deny advisories ===\n' >&2
    nicely cargo deny --locked check advisories "$@" || rc=1
    if command -v bun >/dev/null 2>&1; then
      printf '\n=== bun audit ===\n' >&2
      bun audit || rc=1
    else
      printf 'bun not found; skipping bun audit\n' >&2
    fi
    exit "${rc}"

# Report mise-managed tools that have newer versions available
[group('deps')]
tools-outdated *args:
    #!/usr/bin/env bash
    set -euo pipefail
    mise outdated "$@"

# Regenerate mise.lock from mise.toml (main checkout only)
[group('deps')]
tools-lock *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    require_main_checkout '{{project_dir}}' "tool lockfile updates"
    cd '{{project_dir}}'
    mise install "$@"

# Regenerate lockfiles for the standalone Dylint libraries
[group('deps')]
dylint-lock:
    #!/usr/bin/env bash
    set -euo pipefail
    for manifest in '{{project_dir}}'/dev/lints/{dylint-custom,dylint-mordant,dylint-trailofbits}/Cargo.toml; do
      cargo generate-lockfile --manifest-path "${manifest}"
    done

# Generate the CLI reference docs
[group('docs')]
docs *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    just build "$@"
    '{{project_dir}}/docs/generate-reference-docs' '{{build_dir}}/debug/treetime' '{{project_dir}}/docs/docs/reference.md'
    printf 'Generated: docs/docs/reference.md\n'

# List cargo targets by kind: just list [bin,example,test,bench]
[group('docs')]
list filter="bin,example,test,bench":
    #!/usr/bin/env bash
    set -euo pipefail
    cd '{{project_dir}}'
    cargo metadata --locked --format-version=1 --no-deps 2>/dev/null | python3 -c "
    import json, sys
    data = json.load(sys.stdin)
    kinds = set('{{filter}}'.split(','))
    targets = []
    for pkg in data['packages']:
        for t in pkg['targets']:
            if any(k in kinds for k in t['kind']):
                targets.append((t['kind'][0], t['name']))
    for kind, name in sorted(targets):
        print(f'{kind:<8} {name}')
    "

# ---------------------------------------------------------------------------
# JavaScript / TypeScript (Turbo)
# ---------------------------------------------------------------------------

alias d := desktop-start
alias db := desktop-build
alias a := app-start
alias ab := app-build
alias ar := app-release
alias arw := app-release-watch
alias jl := js-lint
alias jlf := js-lint-fix
alias jf := js-format
alias jfc := js-format-check
alias jc := js-check
alias jt := js-test
alias jot := oxlint-test
alias ji := js-install
alias ju := js-upgrade

# Install JS dependencies (and extract the Electron binary)
[group('js')]
js-install *args:
    #!/usr/bin/env bash
    set -euo pipefail
    cd '{{project_dir}}'
    bun install --frozen-lockfile "$@"
    if [[ ! -f "node_modules/electron/dist/electron" ]]; then
      electron_zip="$(find .cache/electron -maxdepth 1 -name 'electron-v*-linux-x64.zip' -print -quit 2>/dev/null || true)"
      if [[ -n "${electron_zip}" ]]; then
        mkdir -p "node_modules/electron/dist"
        unzip -qo "${electron_zip}" -d "node_modules/electron/dist"
        printf "dist" > "node_modules/electron/path.txt"
      fi
    fi

# Start the desktop dev environment (Vite + Electron + napi)
[group('js')]
desktop-start *args: js-install
    cd '{{project_dir}}' && bun run dev:desktop "$@"

# Build the desktop app
[group('js')]
desktop-build *args: js-install
    cd '{{project_dir}}' && bun run build:desktop "$@"

# Start the web app dev environment (server + Vite)
[group('js')]
app-start *args: js-install
    cd '{{project_dir}}' && bun run dev:web "$@"

# Build the web app
[group('js')]
app-build *args: js-install
    cd '{{project_dir}}' && bun run build:web "$@"

# Release build of the web app (watch)
[group('js')]
app-release *args: js-install
    cd '{{project_dir}}' && bun run release:web "$@"

# Release build of the web app (watch)
[group('js')]
app-release-watch *args: js-install
    cd '{{project_dir}}' && bun run release:web "$@"

# Lint JS/TS (oxlint)
[group('js')]
js-lint *args: js-install
    cd '{{project_dir}}' && bun run lint "$@"

# Lint JS/TS with autofix
[group('js')]
js-lint-fix *args: js-install
    cd '{{project_dir}}' && bun run lint:fix "$@"

# Format JS/TS (oxfmt)
[group('js')]
js-format *args: js-install
    cd '{{project_dir}}' && bun run format "$@"

# Check JS/TS formatting (read-only)
[group('js')]
js-format-check *args: js-install
    cd '{{project_dir}}' && bun run format:check "$@"

# Type-check JS/TS (tsc --noEmit)
[group('js')]
js-check *args: js-install
    cd '{{project_dir}}' && bun run typecheck "$@"

# Run JS/TS tests (vitest)
[group('js')]
js-test *args: js-install
    cd '{{project_dir}}' && bun run test "$@"

# Run custom oxlint rule tests (Node 24; RuleTester is unsupported under Bun)
[group('js')]
oxlint-test *args: js-install
    cd '{{project_dir}}' && node --test "dev/lints/oxlint/__tests__/test_*.ts" "dev/lints/oxlint-anti-slop/**/*.test.ts" "$@"

# Upgrade JS dependencies (main checkout only)
[group('js')]
js-upgrade *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    require_main_checkout '{{project_dir}}' "dependency upgrades"
    cd '{{project_dir}}'
    bun update --latest "$@"

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

# Prepare a fresh checkout or worktree (tools + JS deps)
[group('setup')]
setup:
    #!/usr/bin/env bash
    set -euo pipefail
    cd '{{project_dir}}'
    printf '==> Installing pinned tools (mise)\n'
    if command -v mise >/dev/null 2>&1; then
      mise trust '{{project_dir}}/mise.toml' >/dev/null 2>&1 || true
      mise install
    else
      printf 'mise not found; inside the container tools are prebaked. On a host, install mise: https://mise.jdx.dev\n' >&2
    fi
    printf '==> Installing JS dependencies (bun)\n'
    if command -v bun >/dev/null 2>&1; then
      just js-install
    else
      printf 'bun not found; skipping JS dependency install\n' >&2
    fi
    printf '==> Setup complete\n'

# ---------------------------------------------------------------------------
# App (dev servers)
# ---------------------------------------------------------------------------

# Start the web and API dev servers in the foreground (run under sess)
[group('app')]
up: js-install
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    source '{{project_dir}}/dev/lib/app.sh'
    app_up

# Probe the running dev servers; non-zero when a server is absent or stale
[group('app')]
health:
    #!/usr/bin/env bash
    set -uo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    source '{{project_dir}}/dev/lib/app.sh'
    app_health

# Print resolved ports, the commit each server was built from, and worktree match
[group('app')]
status:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    source '{{project_dir}}/dev/lib/app.sh'
    app_status

# ---------------------------------------------------------------------------
# Generated files
# ---------------------------------------------------------------------------

# Regenerate the OpenAPI document from the Rust server handlers
[group('generated')]
openapi out=(project_dir / "packages/app-contracts/openapi.yaml"):
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{build_dir}}' RUSTFLAGS="$(rustflags_build)"
    kache_use build
    nicely cargo -q run --locked -p app-server --bin generate-openapi -- '{{out}}'

# Regenerate the JSON schemas and fail if any committed copy is stale (read-only)
[group('generated')]
generated-check:
    #!/usr/bin/env bash
    set -uo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    source '{{project_dir}}/dev/lib/generated.sh'
    generated_check

# Report reference/golden fixtures that no test uses, from the registry list
[group('generated')]
fixtures-unused:
    #!/usr/bin/env bash
    set -uo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    source '{{project_dir}}/dev/lib/fixtures.sh'
    fixtures_unused
