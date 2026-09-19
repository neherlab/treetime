#!/usr/bin/env bash
# Report reference and golden fixtures that no test uses (A5).
#
# The authority is the list in dev/registry/reference-files.toml, not a text
# search: every tracked fixture under a registered root must be claimed by an
# entry. This walks the roots, subtracts every listed file, and reports the
# remainder, plus the files the registry itself marks as used by no test.
set -euo pipefail

if ! declare -F project_root >/dev/null 2>&1; then
  # shellcheck source=./utils.sh
  source "${BASH_SOURCE%/*}/utils.sh"
fi

function fixtures_unused() {
  local project_dir registry json rc=0
  project_dir="$(project_root)"
  registry="${project_dir}/dev/registry/reference-files.toml"
  if [[ ! -f "${registry}" ]]; then
    err "fixture registry not found: ${registry}"
    return 1
  fi
  json="$(yq -p toml -o json '.' "${registry}")"

  local roots=() covered=() on_disk=()
  mapfile -t roots < <(jq -r '.roots[]' <<<"${json}")
  mapfile -t covered < <(jq -r '[.reference[].files[]] | unique | .[]' <<<"${json}")
  mapfile -t on_disk < <(git -C "${project_dir}" ls-files -- "${roots[@]}" | grep -v '_capture$' | sort -u)

  # Owning tests and listed files must exist, or the registry is itself stale.
  local name test file
  while IFS= read -r name; do
    test="$(jq -r --arg n "${name}" '.reference[] | select(.name == $n) | .test' <<<"${json}")"
    if [[ -n "${test}" && ! -f "${project_dir}/${test}" ]]; then
      err "group '${name}': owning test not found: ${test}"
      rc=1
    fi
  done < <(jq -r '.reference[].name' <<<"${json}")

  for file in "${covered[@]}"; do
    if [[ ! -f "${project_dir}/${file}" ]]; then
      err "listed fixture is missing on disk: ${file}"
      rc=1
    fi
  done

  # Tracked fixtures that no registry entry claims: candidates for deletion.
  local unlisted=()
  mapfile -t unlisted < <(comm -23 <(printf '%s\n' "${on_disk[@]}") <(printf '%s\n' "${covered[@]}" | sort -u))
  if ((${#unlisted[@]} > 0)); then
    err "tracked fixtures that no registry entry claims (register them in dev/registry/reference-files.toml or delete them):"
    printf '  %s\n' "${unlisted[@]}" >&2
    rc=1
  fi

  # Fixtures the registry records as used by no test, kept on purpose.
  local kept=()
  mapfile -t kept < <(jq -r '.reference[] | select(.used == false) | .files[]' <<<"${json}")
  if ((${#kept[@]} > 0)); then
    info "Fixtures used by no test, kept intentionally (per registry):"
    printf '  %s\n' "${kept[@]}" >&2
  fi

  if ((rc == 0)); then
    info "Every tracked fixture is claimed by the registry."
  fi
  return "${rc}"
}

export -f fixtures_unused
