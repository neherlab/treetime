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

## Proposed implementation plan

The installer does not depend on the open signing decision; the signing steps do.

1. **Installer.** In `packages/app-desktop/electron-builder.config.ts`, set `win.target` to `[{ target: "nsis", arch: ["x64", "arm64"] }]` with `nsis: { oneClick: false, perMachine: false }`, so the app installs per user with a Start menu entry and an uninstaller. Keep `zip` only if a portable download is wanted
2. **Release workflow.** Add a job that runs `dev/desktop/package` for Windows on `windows-2022` and `windows-11-arm` and uploads the installers; no workflow builds the desktop app today
3. **Option A, sign.** With Azure Trusted Signing: set `win.azureSignOptions` in the config and the Azure credentials as CI secrets; sign the CLI `.exe` assets with the same account (`signtool` with the Trusted Signing dlib). With an OV or EV certificate: `win.signtoolOptions` with the certificate from a CI secret or a cloud HSM. Validate with `Get-AuthenticodeSignature` in CI (status `Valid`)
4. **Option B, stay unsigned.** Document the SmartScreen steps (More info, Run anyway) on the installation page
5. **Validate.** Install, start, and uninstall the NSIS installer on both runners in CI (silent mode `/S`), and run one command through the app's N-API addon
