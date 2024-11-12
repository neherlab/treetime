#!/usr/bin/env bash
set -euo pipefail

export JOBS="${JOBS:=$(($(nproc --all) + 2))}"

function print_color() {
  local color_code="${1}"
  local message="${2}"
  if [[ -z "${FORCE_COLOR:-}" && (-n "${NO_COLOR:-}" || "${TERM:-}" == "dumb" || $(tput colors 2>/dev/null) -lt 8) ]]; then
    echo "${message}"
  else
    echo -e "\e[38;5;${color_code}m${message}\e[0m"
  fi
}
export -f print_color

function err() {
  print_color 1 "[ERR] $0: $1"
}
export -f err

function warn() {
  print_color 3 "[WARN] $0: $1"
}
export -f warn

function info() {
  print_color 6 "[INFO] $0: $1"
}
export -f info

function debug() {
  print_color 8 "[DEBUG] $0: $1"
}
export -f debug

function datetime_safe() {
  TZ=UTC date -u '+%Y-%m-%d_%H-%M-%S'
}
export -f datetime_safe

function datetime_ms_safe() {
  TZ=UTC date -u '+%Y-%m-%d_%H-%M-%S_%3NZ'
}
export -f datetime_ms_safe

# Run a command with lower (nicer) CPU and IO priority, time it and print the outcome
function nicely() {
  print_color 8 "+${*:-}"
  if bash -c "
      trap 'print_color 4 \"🟦 Cancelled\" && exit 0' INT;\
      nice -14 ionice -c2 -n3 \
      /usr/bin/time -qf 'Cmd : %C\nTime: %E\nMem : %M KB' \
      ${*}"; then
    (
      { set +x; } 2>/dev/null
      print_color 2 '🟩 Success'
      exit 0
    )
  else
    local cmd_exit_code=$?
    (
      { set +x; } 2>/dev/null
      print_color 1 "🟥 Failure (code: ${cmd_exit_code})"
      exit "${cmd_exit_code}"
    )
  fi
}
export -f nicely

function fake_tty() {
  script -qefc bash -c "$*" /dev/null
}

function file_hash() {
  local paths=("$@")
  local uid=$(id -u)
  local gid=$(id -g)
  local user=$(id -un)
  local group=$(id -gn)
  echo -n "$uid $gid $user $group" | cat - "${paths[@]}" | md5sum | cut -f 1 -d " " | cut -c1-7
}
export -f file_hash

function abspath() {
  readlink -m "${1:?}"
}
export -f abspath

function here() {
  cd -P -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd -P
}
export -f here

function project_root() {
  abspath "$(here)/../.."
}
export -f project_root

function get_build_dir() {
  local target="${1:-}"
  abspath "$(project_root)/.build/docker${target:+"-${target}"}"
}
export -f get_build_dir

function get_cache_dir() {
  local target="${1:-}"
  abspath "$(project_root)/.cache/docker${target:+"-${target}"}"
}
export -f get_cache_dir

function get_test_dir() {
  local target="${1:-}"
  abspath "$(get_build_dir "${target}")/test"
}
export -f get_test_dir

function get_out_dir() {
  abspath "$(project_root)/.out"
}
export -f get_out_dir

function get_bin_dir() {
  local target="${1:-}"
  abspath "$(get_build_dir "${target}")/${target}/release"
}
export -f get_bin_dir

function guess_ext() {
  local target="${1:-}"
  [[ "$target" =~ (mingw|windows) ]] && echo ".exe"
}
export -f guess_ext

function get_full_bin_path() {
  local bin="${1:?}"
  local target="${2:-}"
  ext="$(guess_ext "${target}")"
  echo "$(get_bin_dir "${target}")/${bin}${ext}"
}
export -f get_full_bin_path

function get_final_bin_path() {
  local bin="${1:?}"
  local target="${2:-}"
  ext="$(guess_ext "${target}")"
  echo "$(get_out_dir)/${bin}${target:+-${target}}${ext}"
}
export -f get_final_bin_path

function load_env_maybe() {
  local env_file="${1:-}"
  # shellcheck disable=SC2046
  [ -n "${env_file}" ] && export $(grep -v '^#' .env | xargs)
}
export -f load_env_maybe

function load_env() {
  local env_file="${1:?}"
  if ! [ -r "${env_file}" ]; then
    err "unable to load env file: '${env_file}'" >&2
    exit 1
  fi
  load_env_maybe "${env_file}"
}
export -f load_env

SUCCESS="({ set +x; } 2> /dev/null && echo '🟩 Success' && exit 0)"
FAILURE="({ set +x; } 2> /dev/null && echo '🟥 Failure' && exit 1)"
CANCELLED="({ set +x; } 2> /dev/null && echo '🟦 Cancelled' && exit 0)"
export SUCCESS_OR_FAILURE="&& ${SUCCESS} || ${FAILURE}"
export HANDLE_SIGINT="trap \"${CANCELLED}; exit 0\" INT"
