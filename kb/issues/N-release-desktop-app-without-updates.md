# Desktop app has no update mechanism

## Summary

The desktop app has no update mechanism, so users keep the version they downloaded until they look for a new one.

## Evidence

- [packages/app-desktop/electron-builder.config.ts#L14](../../packages/app-desktop/electron-builder.config.ts#L14): `publish: null`.

> [!IMPORTANT]
> **Decision required.** The update mechanism: `electron-updater` against the GitHub releases (needs signed builds on macOS, [M-release-macos-app-unsigned-and-not-notarized.md](M-release-macos-app-unsigned-and-not-notarized.md)), the package manager of each platform (Flatpak, winget), or a notice in the app that a newer release exists.
