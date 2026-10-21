# macOS app unsigned and not notarized

## Summary

The macOS desktop app is ad-hoc signed, without the hardened runtime and without notarization. Gatekeeper blocks a downloaded app that is not notarized; on macOS 15 the user must allow it in System Settings. The CLI and the app also ship as separate x86_64 and arm64 downloads, so users can pick the wrong one and run it under Rosetta.

## Evidence

- [packages/app-desktop/electron-builder.config.ts#L37-L38](../../packages/app-desktop/electron-builder.config.ts#L37-L38): `identity: "-"`, `hardenedRuntime: false`.
- Release targets `x86_64-apple-darwin` and `aarch64-apple-darwin` in `dev/cross/targets`.

> [!IMPORTANT]
> **Decision required.** Notarization needs an Apple Developer account (yearly fee) with a Developer ID certificate, stored as CI secrets. Options: notarize the app and the CLI; or keep them unsigned and document the steps to open them.

## Potential solutions

- Developer ID signing with the hardened runtime and notarization in the desktop and CLI release workflows.
- Universal2 builds (both architectures in one file) for the CLI and the app.
