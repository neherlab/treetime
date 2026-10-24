#!/usr/bin/env bash
#
# Port slots of the checkouts of one repository. The servers of a checkout use
# the ports of its slot, so the main checkout and its linked worktrees run side
# by side. Portable to bash 3.2 and to the BSD userland of macOS.
#
# The ports lie in the block 14000-14999, in slots of 10 ports: the main checkout
# holds slot 0 (14000-14009), and a linked worktree holds one of the slots 1-99
# (14010-14999). A worktree gets the lowest slot that no other checkout holds and
# whose ports have no listener on 127.0.0.1, and keeps it until the worktree is
# removed. Each claim is a file port-slots/<slot> in the git common directory,
# created exclusively and holding the path of its worktree, so every checkout
# sees the claims of the others and two worktrees never claim one slot. A claim
# of a worktree that no longer exists is deleted before the next slot is chosen.

port_block_start=14000
port_slot_size=10
port_slot_count=100

# Ports in use within a slot, in this order of offsets: api, web, serve, web2,
# web2 serve, docs, docs serve, desktop.
port_slot_used=8

# Print the slot of the checkout that contains the directory $1, and claim a free
# slot when a linked worktree holds none yet.
port_slot() {
  local dir="${1:?}" git_dir common_dir slots_dir toplevel slot
  git_dir="$(git -C "${dir}" rev-parse --path-format=absolute --git-dir)" || return 1
  common_dir="$(git -C "${dir}" rev-parse --path-format=absolute --git-common-dir)" || return 1
  if [[ "${git_dir}" == "${common_dir}" ]]; then
    printf '0\n'
    return 0
  fi
  toplevel="$(cd "$(git -C "${dir}" rev-parse --show-toplevel)" && pwd -P)" || return 1
  slots_dir="${common_dir}/port-slots"

  slot="$(port_slot_find "${slots_dir}" "${toplevel}")"
  if [[ -n "${slot}" ]]; then
    printf '%s\n' "${slot}"
    return 0
  fi

  mkdir -p "${slots_dir}" || return 1
  port_slot_release_removed "${dir}" "${slots_dir}" || return 1
  for ((slot = 1; slot < port_slot_count; slot++)); do
    [[ -e "${slots_dir}/${slot}" ]] && continue
    port_slot_listening "${slot}" && continue
    # noclobber opens the file with O_EXCL, so of two worktrees that race for
    # one slot only the first claims it.
    if (set -o noclobber && printf '%s\n' "${toplevel}" >"${slots_dir}/${slot}") 2>/dev/null; then
      printf '%s\n' "${slot}"
      return 0
    fi
  done
  printf 'No free port slot for %s: all %s worktree slots in %s are taken or listening\n' \
    "${toplevel}" "$((port_slot_count - 1))" "${slots_dir}" >&2
  return 1
}

# Print the port at an offset of a slot.
port_of_slot() {
  printf '%s\n' "$((port_block_start + ${1:?} * port_slot_size + ${2:?}))"
}

# Print the slot that a claim file of the directory $1 assigns to the checkout
# path $2, or nothing.
port_slot_find() {
  local claim
  for claim in "${1:?}"/*; do
    [[ -f "${claim}" ]] || continue
    if [[ "$(cat "${claim}")" == "${2:?}" ]]; then
      printf '%s\n' "${claim##*/}"
      return 0
    fi
  done
}

# Delete the claims in $2 whose worktree no longer exists among the worktrees of
# the repository of the directory $1.
port_slot_release_removed() {
  local dir="${1:?}" slots_dir="${2:?}" worktrees live="" line path claim owner
  worktrees="$(git -C "${dir}" worktree list --porcelain)" || return 1
  while IFS= read -r line; do
    [[ "${line}" == "worktree "* ]] || continue
    path="${line#worktree }"
    if [[ -d "${path}" ]]; then
      live="${live}$(cd "${path}" && pwd -P)"$'\n'
    fi
  done <<<"${worktrees}"
  for claim in "${slots_dir}"/*; do
    [[ -f "${claim}" ]] || continue
    owner="$(cat "${claim}")"
    case $'\n'"${live}" in
    *$'\n'"${owner}"$'\n'*) ;;
    *) rm -f -- "${claim}" ;;
    esac
  done
}

# Succeed when a port of the slot $1 has a listener on 127.0.0.1.
port_slot_listening() {
  local offset port
  for ((offset = 0; offset < port_slot_used; offset++)); do
    port="$(port_of_slot "${1:?}" "${offset}")"
    if (: <>"/dev/tcp/127.0.0.1/${port}") 2>/dev/null; then
      return 0
    fi
  done
  return 1
}
