# TreeTime task runner. Run through the container boundary, e.g.
#   ./dev/docker/run just check
# `just --list` shows every task grouped with a one-line description.
# Tool versions come from mise (mise.toml / mise.lock).

set shell := ["bash", "-euo", "pipefail", "-c"]
set positional-arguments := true

project_dir := justfile_directory()

# Per-kind cargo target directories. The container bind-mounts .build, so these
# match the paths dev/docker/run mounts and keep the build cache warm.
build_dir := project_dir / ".build/docker"
test_dir := build_dir / "test"

# Isolated target dir for the dylint driver: the pinned nightly toolchain must
# not overwrite the stable build cache. The lint library is built under lib/.
dylint_check_dir := build_dir / "dylint/check"
dylint_lib_dir := build_dir / "dylint/lib"

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

# ---------------------------------------------------------------------------
# Lint and format (Rust)
# ---------------------------------------------------------------------------

alias l := lint
alias lf := lint-fix
alias lc := lint-ci
alias lD := lint-deps
alias dy := dylint
alias dyf := dylint-fix
alias f := format
alias fc := format-check
alias q := quality

# Clippy over all targets
[group('lint')]
lint *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{build_dir}}' RUSTFLAGS="$(rustflags_build)"
    kache_use clippy
    nicely cargo -q clippy -q --all-targets --all --locked "$@"

# Clippy with autofix (stage changes first)
[group('lint')]
lint-fix *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    export CARGO_TARGET_DIR='{{build_dir}}' RUSTFLAGS="$(rustflags_build)"
    kache_use clippy
    vcs_flag="--allow-staged"; [[ -f '{{project_dir}}/.git' ]] && vcs_flag="--allow-no-vcs"
    nicely cargo -q clippy -q --all-targets --all --fix "${vcs_flag}" --locked "$@"

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

# Run the custom dylint lint library
[group('lint')]
dylint *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    so="$(dylint_build_lib '{{project_dir}}' '{{dylint_lib_dir}}')"
    export CARGO_TARGET_DIR='{{dylint_check_dir}}' RUSTFLAGS="$(rustflags_build)" RUST_BACKTRACE=0
    kache_use dylint
    nicely cargo dylint --quiet --lib-path "${so}" -- --quiet --locked --workspace --all-targets "$@"

# Run the custom dylint lint library with autofix
[group('lint')]
dylint-fix *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    so="$(dylint_build_lib '{{project_dir}}' '{{dylint_lib_dir}}')"
    export CARGO_TARGET_DIR='{{dylint_check_dir}}' RUSTFLAGS="$(rustflags_build)" RUST_BACKTRACE=0
    kache_use dylint
    vcs_flag="--allow-staged"; [[ -f '{{project_dir}}/.git' ]] && vcs_flag="--allow-no-vcs"
    nicely cargo dylint --quiet --fix --lib-path "${so}" -- "${vcs_flag}" --quiet --locked --workspace --all-targets "$@"

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
check-full: (_check "full")

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
    run_check "rust-clippy" cargo -q clippy -q --all-targets --all --locked
    if have_bun_project; then
      run_check "typescript" bash -c "cd '{{project_dir}}' && bun run typecheck"
      run_check "oxlint" bash -c "cd '{{project_dir}}' && bun run lint"
    else
      skip "typescript" "no node_modules; run just setup"
      skip "oxlint" "no node_modules; run just setup"
    fi

    if [[ "{{mode}}" == "full" ]]; then
      if command -v shellcheck >/dev/null 2>&1; then
        mapfile -t sh_files < <(dev_shell_files '{{project_dir}}')
        run_check "shellcheck" shellcheck "${sh_files[@]}"
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
      if command -v cargo-dylint >/dev/null 2>&1; then
        run_check "dylint" bash -c "just dylint"
      else
        skip "dylint" "cargo-dylint not installed"
      fi
      if [[ -f '{{project_dir}}/deny.toml' ]] && command -v cargo-deny >/dev/null 2>&1; then
        run_check "cargo-deny" cargo deny --locked check
      else
        skip "cargo-deny" "no deny.toml"
      fi
      if command -v cargo-shear >/dev/null 2>&1; then
        run_check "unused-deps" cargo shear
      else
        skip "unused-deps" "cargo-shear not installed"
      fi
      if have_bun_project && grep -q '"knip"' package.json 2>/dev/null; then
        run_check "knip" bash -c "cd '{{project_dir}}' && bun run knip"
      else
        skip "knip" "not configured"
      fi
      skip "generated-freshness" "no freshness command yet"
      export CARGO_TARGET_DIR='{{test_dir}}' RUSTFLAGS="$(rustflags_test)"
      export KACHE_CACHE_DIR="${kache_base}"; kache_use test
      run_check "tests" cargo -q nextest run --locked --workspace --cargo-quiet --no-fail-fast --hide-progress-bar
      if have_bun_project; then
        run_check "js-tests" bash -c "cd '{{project_dir}}' && bun run test"
      else
        skip "js-tests" "no node_modules"
      fi
    fi

    printf '\n===== %s check summary =====\n' "{{mode}}"
    fail=0
    for i in "${!names[@]}"; do
      printf '  %-22s %s\n' "${names[$i]}" "${results[$i]}"
      [[ "${results[$i]}" == FAIL* ]] && fail=$((fail + 1))
    done
    printf '\n'
    if (( fail > 0 )); then
      printf 'Reported %d failing check(s). This gate is report-only and does not fail the command.\n' "${fail}"
    else
      printf 'All checks that ran passed.\n'
    fi

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
    nicely shellcheck "$@" "${files[@]}"

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

# Explain why a crate is in the dependency tree: just why <crate>
[group('deps')]
why *args:
    #!/usr/bin/env bash
    set -euo pipefail
    source '{{project_dir}}/dev/lib/utils.sh'
    nicely cargo -q tree -i -p --locked "$@"

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
alias ji := js-install
alias ju := js-upgrade

# Install JS dependencies (and extract the Electron binary)
[group('js')]
js-install *args:
    #!/usr/bin/env bash
    set -euo pipefail
    cd '{{project_dir}}'
    bun install "$@"
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
