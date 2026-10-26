# Windows app unsigned and without installer

## Summary

The Windows desktop app is packaged as a zip archive without Authenticode signing. Users get no Start menu entry, no uninstaller, and a SmartScreen warning about an unknown publisher.

## Evidence

- [packages/app-desktop/electron-builder.config.ts#L32](../../packages/app-desktop/electron-builder.config.ts#L32): `target: "zip"`.
- No signing configuration in the desktop packaging (`dev/desktop/package`, `packages/app-desktop/electron-builder.config.ts`), and no workflow builds or publishes the desktop app.

> [!IMPORTANT]
> **Decision required.** Code signing needs a certificate: an OV or EV certificate, or Azure Trusted Signing (monthly fee). Options: sign the installer and the CLI; or keep them unsigned and document the SmartScreen steps.

## Potential solutions

- An NSIS or MSIX installer from electron-builder, signed, for x64 and arm64.
