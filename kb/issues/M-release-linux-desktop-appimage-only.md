# Linux desktop app ships as AppImage only

## Summary

The Linux desktop app ships only as an AppImage. Electron AppImages are reported to fail on Ubuntu 24.04 and later, where AppArmor restricts the unprivileged user namespaces that the Chromium sandbox needs, and older AppImage runtimes need `libfuse2`, which current distributions do not install by default.

## Evidence

- [packages/app-desktop/electron-builder.config.ts#L23-L27](../../packages/app-desktop/electron-builder.config.ts#L23-L27): AppImage toolset 1.0.3, `target: "AppImage"`.

> [!IMPORTANT]
> **Investigation required.** Start the AppImage of the newest nightly (toolset 1.0.3) on a fresh Ubuntu 24.04 desktop and on a fresh Fedora desktop, both without `libfuse2`. Record the exit message, whether the Chromium sandbox or FUSE stops it, and whether `--no-sandbox` or an AppArmor profile makes it start.

> [!IMPORTANT]
> **Decision required.** The Linux package formats: Flatpak on Flathub (sandbox, updates, discovery), `.deb` and `.rpm` from electron-builder, or both, next to or instead of the AppImage.
