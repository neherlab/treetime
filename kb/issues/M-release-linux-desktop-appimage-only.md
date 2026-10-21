# Linux desktop app ships as AppImage only

## Summary

The Linux desktop app ships only as an AppImage. Electron AppImages are reported to fail on Ubuntu 24.04 and later, where AppArmor restricts the unprivileged user namespaces that the Chromium sandbox needs, and older AppImage runtimes need `libfuse2`, which current distributions do not install by default.

## Evidence

- [packages/app-desktop/electron-builder.config.ts#L23-L27](../../packages/app-desktop/electron-builder.config.ts#L23-L27): AppImage toolset 1.0.3, `target: "AppImage"`.

> [!IMPORTANT]
> **Investigation required.** Start the AppImage on a fresh Ubuntu 24.04 desktop and on Fedora, and record whether the sandbox or FUSE stops it with toolset 1.0.3.

> [!IMPORTANT]
> **Decision required.** The Linux package formats: Flatpak on Flathub (sandbox, updates, discovery), `.deb` and `.rpm` from electron-builder, or both, next to or instead of the AppImage.
