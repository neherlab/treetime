# JavaScript dev shortcuts fail without the Electron cache

## Summary

Every JavaScript shortcut routed through `bun_run()` invokes `js-install` before its requested command. When `node_modules/electron/dist/electron` and `.cache/electron` are both absent, `js-install` exits with status 1 while searching for a cached Electron archive. The requested command never starts, including commands that do not use Electron such as Oxlint and Oxformat.

## Reproduction

From a clean worktree without `.cache/electron`:

```bash
./dev/docker/run ./dev/dev jl
```

`bun install` completes, then the wrapper exits with status 1 before printing or executing `bun run lint`. Repeating the command after dependencies are installed produces the same exit status.

## Cause

`bun_run()` unconditionally calls `run_one js-install` before dispatching any JavaScript package script [`dev/dev#L92-L97`](../../dev/dev#L92-L97). When the Electron binary is absent, `js-install` executes `find .cache/electron ...` in a command substitution under `set -e` [`dev/dev#L279-L293`](../../dev/dev#L279-L293). `find` returns a non-zero status when the cache directory does not exist, so the shell exits before the optional archive check or requested JavaScript command.

## Impact

Fresh worktrees cannot run JavaScript linting, formatting, type checks, tests, web builds, or desktop builds through the required `./dev/docker/run ./dev/dev` interface unless the Electron cache directory already exists.

## Potential solutions

- Guard the archive search with a directory-existence check and treat an absent optional cache as an empty result.
- Move Electron binary preparation to desktop commands so non-Electron JavaScript commands do not depend on desktop runtime state.

The first option fixes the immediate shell failure. The second also separates general JavaScript tooling from the desktop runtime contract; desktop commands would still need an explicit error when neither an installed binary nor cached archive is available.
