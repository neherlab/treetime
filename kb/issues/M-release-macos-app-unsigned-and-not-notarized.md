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

## Proposed implementation plan

No universal2 build exists today (no `lipo` or `universal` target in `dev/`, `.github/`, or the desktop config). Signing depends on the open decision above; the universal2 build does not.

1. **Universal2 builds.** CLI: after the two darwin builds in the release workflow, join them with `lipo -create` into one asset and upload it instead of the two. App: set `mac.target` in `packages/app-desktop/electron-builder.config.ts` to `{ target: "dmg", arch: ["universal"] }`; the `app-napi` addon then needs both architectures built and merged (`@electron/universal` merges them when both `.node` files are present). Validate with `lipo -archs` on the outputs and a run on `macos-15-intel` and an Apple Silicon runner in `compat-macos`
2. **Option A, notarize.** Store the Developer ID certificate and an App Store Connect API key as CI secrets; set `identity` to the Developer ID name, `hardenedRuntime: true`, an entitlements file with `com.apple.security.cs.allow-jit` for V8, and `notarize: true`. Sign the CLI with `codesign --options runtime` and notarize it with `xcrun notarytool submit --wait` (zipped). Validate with `spctl --assess --type execute` and `xcrun stapler validate` in CI
3. **Option B, stay unsigned.** Document on the installation page how to open the app (System Settings, Privacy and Security, Open Anyway) and how to clear the quarantine attribute of the CLI (`xattr -d com.apple.quarantine`)
