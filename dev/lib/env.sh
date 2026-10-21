#!/usr/bin/env bash
#
# Shell helpers shared by the scripts in dev/. Portable to bash 3.2 and to the
# BSD userland of macOS.

# Export KEY=VALUE pairs of a dotenv file. A variable that is already set in the
# environment wins over the file. Blank lines, comment lines, an `export `
# prefix, surrounding quotes, and a trailing ` # comment` after an unquoted value
# are handled; the file is parsed, never executed.
load_env_file() {
  local file="${1:?}" line key value
  [[ -r "${file}" ]] || return 0
  while IFS= read -r line || [[ -n "${line}" ]]; do
    line="${line#"${line%%[![:space:]]*}"}"
    [[ -z "${line}" || "${line}" == \#* ]] && continue
    line="${line#export }"
    [[ "${line}" == *=* ]] || continue
    key="${line%%=*}"
    value="${line#*=}"
    [[ "${key}" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] || continue
    case "${value}" in
    \"*\"*) value="${value#\"}" && value="${value%%\"*}" ;;
    \'*\'*) value="${value#\'}" && value="${value%%\'*}" ;;
    *) value="${value%%[[:space:]]#*}" && value="${value%"${value##*[![:space:]]}"}" ;;
    esac
    if [[ -z "${!key+x}" ]]; then
      export "${key}=${value}"
    fi
  done <"${file}"
}

# Print a short content hash of files and directories (recursively), their
# relative paths included, plus any extra literal strings given after them.
hash_files() {
  local sha
  if command -v sha256sum >/dev/null 2>&1; then
    sha=(sha256sum)
  else
    sha=(shasum -a 256)
  fi
  local arg
  {
    for arg in "$@"; do
      if [[ -d "${arg}" ]]; then
        find "${arg}" -type f | LC_ALL=C sort | while IFS= read -r file; do
          "${sha[@]}" "${file}"
        done
      elif [[ -f "${arg}" ]]; then
        "${sha[@]}" "${arg}"
      else
        printf '%s\n' "${arg}"
      fi
    done
  } | "${sha[@]}" | cut -c1-16
}

# Print the CPU budget of the containers of dev/docker/run: TREETIME_DOCKER_CPUS,
# default 8, lowered to the CPU count of the Docker daemon, because `docker run
# --cpus` rejects a value above it.
docker_cpu_budget() {
  local budget="${TREETIME_DOCKER_CPUS:-8}" available
  if ! [[ "${budget}" =~ ^[0-9]+([.][0-9]+)?$ ]] || ! awk -v b="${budget}" 'BEGIN { exit !(b > 0) }'; then
    printf 'TREETIME_DOCKER_CPUS must be a positive number of CPUs, got: %s\n' "${budget}" >&2
    return 1
  fi
  available="$(docker info --format '{{.NCPU}}' 2>/dev/null)" || available=""
  awk -v b="${budget}" -v a="${available}" 'BEGIN { if (a > 0 && b > a) b = a; print b }'
}

# Print the share of a CPU budget for one of several containers that run at the
# same time, with two decimals and at least the 0.01 CPUs that Docker accepts.
docker_cpu_share() {
  awk -v b="${1:?}" -v n="${2:?}" 'BEGIN { c = b / n; if (c < 0.01) c = 0.01; printf "%.2f\n", c }'
}
