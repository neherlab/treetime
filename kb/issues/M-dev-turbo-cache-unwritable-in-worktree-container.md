# Turbo cache is unwritable in the container of a linked worktree

`just typecheck` fails in the build container of a linked git worktree, so `just check` and `just check-all` fail there on a clean tree:

```
x IO error: failed to create directory `<main-checkout>/.turbo/cache`
`-> Permission denied (os error 13)
```

## Cause

Turborepo shares its local cache with the main worktree by default: in a linked worktree, `turbo run` writes to `.turbo/cache` of the main checkout ([turborepo.dev: cacheDir](https://turborepo.dev/docs/reference/configuration)). `dev/docker/run` mounts the current checkout read-write and the git common directory read-only, and the main checkout's `.turbo/` is not writable from the container of a linked worktree.

`turbo.json` sets no `cacheDir`, so the default sharing applies. On the host the shared directory is writable and the failure does not occur.

The turbo documentation also notes that a shared cache hit can restore outputs that contain absolute paths of another worktree.

## Workaround

Set a worktree-local cache for the command:

```bash
./dev/docker/run bash -c 'TURBO_CACHE_DIR=.turbo/cache just check-all'
```

## Open question

- **Per-worktree cache**: set `"cacheDir": ".turbo/cache"` in `turbo.json`. A relative `cacheDir` resolves from the current worktree and disables sharing, so each worktree builds its own cache, on the host too
- **Shared cache**: mount the main checkout's `.turbo/` read-write in `dev/docker/run` when it runs in a linked worktree, keeping cross-worktree cache hits and the absolute-path caveat above
