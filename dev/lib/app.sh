#!/usr/bin/env bash
# App dev-server orchestration for the `up`, `health`, and `status` tasks.
#
# `up` runs the API server and the Vite dev server together in the foreground so
# it lives under `sess` like any long-running command; `health` and `status` are
# read-only probes of what `up` recorded. Ports are resolved per worktree so two
# checkouts can run side by side without colliding, and the real environment wins
# over the resolved value and over `.env.local`.
set -euo pipefail

if ! declare -F project_root >/dev/null 2>&1; then
  # shellcheck source=./utils.sh
  source "${BASH_SOURCE%/*}/utils.sh"
fi

# Base ports. The API server listens here and Vite proxies /api to it; Vite
# serves the browser app on the web port.
APP_API_BASE_PORT=3100
APP_WEB_BASE_PORT=5173

APP_API_PID=""
APP_WEB_PID=""

# A deterministic 0..511 offset unique to this worktree, derived from its root
# path so each checkout resolves to its own port pair.
function app_port_offset() {
  local root sum
  root="$(git -C "$(project_root)" rev-parse --show-toplevel 2>/dev/null || project_root)"
  sum="$(printf '%s' "${root}" | cksum | cut -d' ' -f1)"
  printf '%s' "$((sum % 512))"
}

# Resolve API_PORT and WEB_PORT for this worktree. An explicit environment value
# wins over the resolved one: PORT or TREETIME_API_PORT for the API, and
# TREETIME_WEB_PORT for the web server.
function app_resolve_ports() {
  local offset
  offset="$(app_port_offset)"
  API_PORT="${TREETIME_API_PORT:-${PORT:-$((APP_API_BASE_PORT + offset))}}"
  WEB_PORT="${TREETIME_WEB_PORT:-$((APP_WEB_BASE_PORT + offset))}"
  export API_PORT WEB_PORT
}

function app_commit() {
  git -C "$(project_root)" rev-parse --short HEAD 2>/dev/null || printf 'unknown'
}

function app_dirty() {
  if [[ -n "$(git -C "$(project_root)" status --porcelain 2>/dev/null)" ]]; then
    printf '1'
  else
    printf '0'
  fi
}

function app_state_dir() {
  printf '%s/tmp/up' "$(project_root)"
}

function app_state_file() {
  printf '%s/%s.json' "$(app_state_dir)" "${1:?service name required}"
}

# Prefix every line read on stdin with a service label, so both servers can share
# one terminal legibly.
function app_prefix() {
  local label="${1:?label required}" line
  while IFS= read -r line; do
    printf '[%s] %s\n' "${label}" "${line}"
  done
}

# GET a URL, discarding the body; non-zero when the server does not answer.
function app_probe() {
  curl -fsS -o /dev/null --max-time 3 "${1:?url required}"
}

# Poll a health URL until it answers or the timeout elapses. Probing the socket
# is the only signal a separate process exposes for "ready"; the timeout is a
# safety bound, not a fixed synchronization delay.
function app_wait_ready() {
  local url="${1:?url required}" timeout="${2:-180}" elapsed=0
  while ((elapsed < timeout)); do
    if app_probe "${url}"; then
      return 0
    fi
    sleep 1
    elapsed=$((elapsed + 1))
  done
  return 1
}

# Record what a server was started as, so health and status can compare it later.
function app_write_state() {
  local svc="${1:?}" url="${2:?}" port="${3:?}"
  mkdir -p "$(app_state_dir)"
  jq -n \
    --arg svc "${svc}" \
    --arg url "${url}" \
    --argjson port "${port}" \
    --arg commit "$(app_commit)" \
    --argjson dirty "$(app_dirty)" \
    --arg started "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" \
    '{service: $svc, url: $url, port: $port, commit: $commit, dirty: ($dirty == 1), started_at: $started}' \
    >"$(app_state_file "${svc}")"
}

function app_cleanup() {
  trap - INT TERM EXIT
  local pid
  for pid in "${APP_WEB_PID}" "${APP_API_PID}"; do
    [[ -n "${pid}" ]] && kill -TERM -"${pid}" 2>/dev/null || true
  done
}

# Start both servers in order, stream their output, print one ready line, and
# block until Ctrl-C tears the whole tree down.
function app_up() {
  local project_dir
  project_dir="$(project_root)"
  load_env_maybe "${project_dir}/.env.local"
  app_resolve_ports
  mkdir -p "$(app_state_dir)"

  local api_url="http://127.0.0.1:${API_PORT}/api/version"
  local web_url="http://127.0.0.1:${WEB_PORT}/"
  local api_cmd web_cmd
  api_cmd="cd '${project_dir}' && PORT='${API_PORT}' exec just r treetime-server -- --data-dir=data --out-dir=tmp/web"
  web_cmd="cd '${project_dir}/packages/app-web' && PORT='${API_PORT}' exec bun run dev:nowatch -- --port '${WEB_PORT}'"

  trap app_cleanup INT TERM EXIT

  info "Starting API server on http://127.0.0.1:${API_PORT} (web proxies /api here)"
  setsid bash -c "${api_cmd}" > >(app_prefix api) 2>&1 &
  APP_API_PID=$!
  # The first probe waits through a possible cold compile of the server binary.
  if ! app_wait_ready "${api_url}" 900; then
    err "API server did not answer at ${api_url}"
    return 1
  fi
  app_write_state api "${api_url}" "${API_PORT}"

  info "Starting web server on http://127.0.0.1:${WEB_PORT}"
  setsid bash -c "${web_cmd}" > >(app_prefix web) 2>&1 &
  APP_WEB_PID=$!
  if ! app_wait_ready "${web_url}" 300; then
    err "Web server did not answer at ${web_url}"
    return 1
  fi
  app_write_state web "${web_url}" "${WEB_PORT}"

  print_color 2 "Ready: web http://127.0.0.1:${WEB_PORT}  |  api http://127.0.0.1:${API_PORT}/api  (Ctrl-C to stop)"
  wait "${APP_API_PID}" "${APP_WEB_PID}"
}

# Probe both servers. Non-zero when a server is absent (never started or not
# answering) or stale (built from a commit other than the current worktree HEAD).
function app_health() {
  local rc=0 svc state url built commit_now
  commit_now="$(app_commit)"
  for svc in api web; do
    state="$(app_state_file "${svc}")"
    if [[ ! -f "${state}" ]]; then
      err "${svc}: not running (no state; start it with 'just up')"
      rc=1
      continue
    fi
    url="$(jq -r '.url' "${state}")"
    if ! app_probe "${url}"; then
      err "${svc}: not answering at ${url}"
      rc=1
      continue
    fi
    built="$(jq -r '.commit' "${state}")"
    if [[ "${built}" != "${commit_now}" ]]; then
      err "${svc}: stale (built from ${built}, worktree at ${commit_now})"
      rc=1
      continue
    fi
    info "${svc}: ok (${url}, built from ${built})"
  done
  return "${rc}"
}

# Print the resolved ports, the commit each server was started from, and whether
# that matches the current worktree.
function app_status() {
  app_resolve_ports
  local commit_now dirty
  commit_now="$(app_commit)"
  dirty="$(app_dirty)"

  printf 'Resolved ports (this worktree):\n'
  printf '  api  %s\n  web  %s\n' "${API_PORT}" "${WEB_PORT}"
  printf 'Worktree HEAD: %s%s\n' "${commit_now}" "$([[ "${dirty}" == 1 ]] && printf ' (dirty)')"
  printf 'Servers:\n'

  local svc state built surl match
  for svc in api web; do
    state="$(app_state_file "${svc}")"
    if [[ ! -f "${state}" ]]; then
      printf '  %-4s not started (no state)\n' "${svc}"
      continue
    fi
    built="$(jq -r '.commit' "${state}")"
    surl="$(jq -r '.url' "${state}")"
    if [[ "${built}" == "${commit_now}" && "${dirty}" == 0 ]]; then
      match="yes"
    else
      match="no"
    fi
    printf '  %-4s built from %s  url %s  matches worktree: %s\n' "${svc}" "${built}" "${surl}" "${match}"
  done
}

export -f app_port_offset app_resolve_ports app_commit app_dirty
export -f app_state_dir app_state_file app_prefix app_probe app_wait_ready
export -f app_write_state app_cleanup app_up app_health app_status
