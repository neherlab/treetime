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

# Run a command, propagate signals to entire process tree.
# setsid: own process group so group-kill reaches all descendants.
# bg+wait: so bash can deliver trap handlers during long-running children.
# Double-wait: retrieves real exit code after signal interruption.
# See: https://veithen.io/2014/11/16/sigterm-propagation.html
function run_with_signals() {
  setsid "$@" &
  local pid=$!
  trap 'kill -TERM -${pid} 2>/dev/null' INT TERM
  wait "${pid}"
  trap - INT TERM
  wait "${pid}"
}
export -f run_with_signals

# Run a command with lower (nicer) CPU and IO priority, time it and print the outcome
function nicely() {
  local -a env_vars=()
  while (( $# > 0 )) && [[ "${1}" == *=* ]]; do
    env_vars+=("${1}")
    shift
  done

  print_color 8 "+${*:-}"

  local key
  local -a saved_keys=()
  local -A saved_vals=()
  local -A saved_set=()
  for pair in "${env_vars[@]}"; do
    key="${pair%%=*}"
    saved_keys+=("${key}")
    if [[ -v "${key}" ]]; then
      saved_set["${key}"]=1
      saved_vals["${key}"]="${!key}"
    fi
    export "${pair}"
  done

  local cmd_exit_code=0
  run_with_signals \
    nice -14 ionice -c2 -n3 \
    /usr/bin/time -qf 'Cmd : %C\nTime: %E\nMem : %M KB' \
    "$@" || cmd_exit_code=$?

  for key in "${saved_keys[@]}"; do
    if [[ "${saved_set["${key}"]:-}" == "1" ]]; then
      export "${key}=${saved_vals["${key}"]}"
    else
      unset "${key}"
    fi
  done
  if [[ ${cmd_exit_code} -eq 0 ]]; then
    print_color 2 '🟩 Success'
  elif [[ ${cmd_exit_code} -gt 128 && ${cmd_exit_code} -le 159 ]]; then
    # 128+N where N is a signal (1-31): killed by signal
    print_color 4 "🟦 Cancelled"
  else
    print_color 1 "🟥 Failure (code: ${cmd_exit_code})"
    return "${cmd_exit_code}"
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

function copy_bin_to_out() {
  local bin="${1:?}"
  local target="${2:-}"
  local full_bin="$(get_full_bin_path "${bin}" "${target}")"
  local final_bin="$(get_final_bin_path "${bin}" "${target}")"
  mkdir -p "$(dirname "${final_bin}")"
  cp "${full_bin}" "${final_bin}"
}
export -f copy_bin_to_out

function load_env_maybe() {
  local env_file="${1:-}"
  [[ -n "${env_file}" && -r "${env_file}" ]] || return 0
  local line key
  while IFS= read -r line || [[ -n "${line}" ]]; do
    key="${line%%=*}"
    [[ -z "${!key+x}" ]] && export "${line}"
  done < <(grep -vE '^(#|$)' "${env_file}") || true
}
export -f load_env_maybe

function load_env() {
  local env_file="${1:?}"
  if ! [[ -r "${env_file}" ]]; then
    err "unable to load env file: '${env_file}'" >&2
    exit 1
  fi
  load_env_maybe "${env_file}"
}
export -f load_env

function getenv_cross() {
  local varname="${1:?}"
  local target="${2:-${CROSS_COMPILE:?}}"
  printenv "${varname^^}_${target}"
}
export -f getenv_cross

SUCCESS="({ set +x; } 2> /dev/null && echo '🟩 Success' && exit 0)"
FAILURE="({ set +x; } 2> /dev/null && echo '🟥 Failure' && exit 1)"
CANCELLED="({ set +x; } 2> /dev/null && echo '🟦 Cancelled' && exit 0)"
export SUCCESS_OR_FAILURE="&& ${SUCCESS} || ${FAILURE}"
export HANDLE_SIGINT="trap \"${CANCELLED}; exit 0\" INT"
