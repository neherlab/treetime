# Linux desktop app ships as AppImage only

## Summary

The Linux desktop app ships only as an AppImage. Electron AppImages are reported to fail on Ubuntu 24.04 and later, where AppArmor restricts the unprivileged user namespaces that the Chromium sandbox needs, and older AppImage runtimes need `libfuse2`, which current distributions do not install by default.

## Evidence

- [packages/app-desktop/electron-builder.config.ts#L23-L27](../../packages/app-desktop/electron-builder.config.ts#L23-L27): AppImage toolset 1.0.3, `target: "AppImage"`.

> [!IMPORTANT]
> **Investigation required.** Start the AppImage on a fresh Ubuntu 24.04 desktop and on Fedora, and record whether the sandbox or FUSE stops it with toolset 1.0.3.

> [!IMPORTANT]
> **Decision required.** The Linux package formats: Flatpak on Flathub (sandbox, updates, discovery), `.deb` and `.rpm` from electron-builder, or both, next to or instead of the AppImage.

## Proposed implementation plan

The plan depends on the investigation and on the open decision above.

1. **Investigate.** Start the AppImage of the newest nightly on a fresh Ubuntu 24.04 desktop VM and a fresh Fedora VM, without `libfuse2`; record the exit message and whether `--no-sandbox` or an AppArmor profile is needed
2. **Option A, `.deb` and `.rpm`**: add `deb` and `rpm` to `linux.target` in `packages/app-desktop/electron-builder.config.ts`; the `.deb` installs an AppArmor profile that allows user namespaces for the app binary (as Chrome and VS Code do). Upload both from the release workflow. Validate by installing each package in the VMs of step 1 and running one command
3. **Option B, Flatpak on Flathub**: add a Flatpak manifest based on `org.electronjs.Electron2.BaseApp`, build it in CI with `flatpak-builder`, and submit it to Flathub. Validate with `flatpak run` in the VMs of step 1
4. **Docs.** Update the desktop section of the installation page (see [M-release-install-docs-and-system-requirements-missing.md](M-release-install-docs-and-system-requirements-missing.md)) with the chosen formats
