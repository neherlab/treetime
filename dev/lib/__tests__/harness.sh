#!/usr/bin/env bash
#
# Test runner of the shell libraries in dev/lib, sourced by the test scripts in
# this directory. A test is a function whose name starts with the prefix that
# the script passes to run_tests. Each test runs in its own subshell and
# temporary directory, with the library under test sourced.

# Run every test function whose name starts with a prefix, with a library
# sourced, and print one line per test. Return 1 when a test fails.
run_tests() {
  local prefix="${1:?prefix}" library="${2:?library}" name status=0
  for name in $(compgen -A function "${prefix}"); do
    if run_test "${name}" "${library}"; then
      printf 'ok    %s\n' "${name}"
    else
      printf 'FAIL  %s\n' "${name}"
      status=1
    fi
  done
  return "${status}"
}

run_test() {
  local name="${1}" library="${2}" work_dir status=0
  work_dir="$(mktemp -d)"
  (
    cd "${work_dir}"
    # shellcheck source=/dev/null
    source "${library}"
    "${name}"
  ) || status=$?
  rm -rf "${work_dir}"
  return "${status}"
}

assert_eq() {
  [[ "${1}" == "${2}" ]] || fail "expected '${1}', got '${2}'"
}

fail() {
  printf '  %s\n' "${1}" >&2
  exit 1
}
