#!/usr/bin/env bash
# Freshness check for generated files.
#
# The JSON schemas under packages/schemas are produced by `treetime schema`.
# This regenerates them into a temporary directory and compares byte for byte
# against the committed copies, failing on any drift so a stale schema cannot be
# committed. Nothing here writes into the working tree.
set -euo pipefail

if ! declare -F project_root >/dev/null 2>&1; then
  # shellcheck source=./utils.sh
  source "${BASH_SOURCE%/*}/utils.sh"
fi

function generated_check() {
  local project_dir committed tmp rc=0
  project_dir="$(project_root)"
  committed="${project_dir}/packages/schemas"
  tmp="$(mktemp -d)"
  # shellcheck disable=SC2064
  trap "rm -rf '${tmp}'" RETURN

  info "Regenerating JSON schemas to compare against committed copies"
  if ! (cd "${project_dir}" && just r treetime -- schema --for all -o "${tmp}") >&2; then
    err "schema generation failed"
    return 1
  fi

  local file name
  for file in "${committed}"/*.schema.json; do
    name="$(basename "${file}")"
    if [[ ! -f "${tmp}/${name}" ]]; then
      err "committed schema ${name} is no longer produced by the generator"
      rc=1
      continue
    fi
    if ! diff -u "${file}" "${tmp}/${name}" >&2; then
      err "committed schema ${name} is stale; regenerate with 'just r treetime -- schema --for all -o packages/schemas'"
      rc=1
    fi
  done

  for file in "${tmp}"/*.schema.json; do
    name="$(basename "${file}")"
    if [[ ! -f "${committed}/${name}" ]]; then
      err "generator produced ${name}, which has no committed copy in packages/schemas"
      rc=1
    fi
  done

  local openapi_committed openapi_tmp
  openapi_committed="${project_dir}/packages/app-contracts/openapi.yaml"
  openapi_tmp="${tmp}/openapi.yaml"
  info "Regenerating the OpenAPI document to compare against the committed copy"
  if ! (cd "${project_dir}" && just openapi "${openapi_tmp}") >&2; then
    err "OpenAPI generation failed"
    rc=1
  elif ! diff -u "${openapi_committed}" "${openapi_tmp}" >&2; then
    err "committed openapi.yaml is stale; regenerate with 'just openapi'"
    rc=1
  fi

  if ((rc == 0)); then
    info "All generated schema files are up to date."
  fi
  return "${rc}"
}

export -f generated_check
