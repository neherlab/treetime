# Desktop app has no update mechanism

## Summary

The desktop app has no update mechanism, so users keep the version they downloaded until they look for a new one.

## Evidence

- [packages/app-desktop/electron-builder.config.ts#L14](../../packages/app-desktop/electron-builder.config.ts#L14): `publish: null`.

> [!IMPORTANT]
> **Decision required.** The update mechanism: `electron-updater` against the GitHub releases (needs signed builds on macOS, [M-release-macos-app-unsigned-and-not-notarized.md](M-release-macos-app-unsigned-and-not-notarized.md)), the package manager of each platform (Flatpak, winget), or a notice in the app that a newer release exists.

## Proposed implementation plan

No code uses `electron-updater` or another update library in `packages/app-desktop`. The plan depends on the open decision above; one short plan per option:

- **Option A, `electron-updater` on GitHub releases**: set `publish` in `packages/app-desktop/electron-builder.config.ts` to the GitHub provider, upload the `latest*.yml` files with the release assets in the release workflow, and call `autoUpdater.checkForUpdatesAndNotify()` in the main process after the window opens. Needs signed macOS builds first. Test: a unit test of the pure function that decides whether to offer an update from the current and the feed version; validate by hand with a local feed of a higher version
- **Option B, platform package managers**: no app code; depends on [M-release-package-channels-missing.md](M-release-package-channels-missing.md)
- **Option C, notice in the app**: on start, fetch the latest release tag from the GitHub API in the main process, compare it with the app version (`semver`), and show a dismissible notice with a link to the release page. Test: unit test of the version comparison, with release tags of the repository as oracle; MSW mock of the API for the notice component

Validate with `just typecheck`, `just lint-ts`, `just test-ts`.
