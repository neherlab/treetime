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
