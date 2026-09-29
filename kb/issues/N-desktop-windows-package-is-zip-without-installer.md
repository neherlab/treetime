# The Windows desktop package is a zip archive without an installer

The Windows package of the desktop app (`treetime-desktop-x86_64-pc-windows-gnu.zip`) is a zip archive of the unpacked app. Users extract it and start `TreeTime.exe`; nothing adds a Start menu entry or an uninstaller. An NSIS installer (electron-builder target `nsis`) would install the app for the current user and register it with Windows.

## Evidence

- electron-builder compiles the NSIS installer with its own `makensis`, then runs the finished installer once under Wine to create the uninstaller (`computeScriptAndSignUninstaller` in `packages/app-builder-lib/src/targets/nsis/NsisTarget.ts`)
- On Linux, electron-builder 26.16.1 with `toolsets.wine: "1.0.1"` downloads `wine-11.0-linux-x86_64.tar.xz` from the release `wine@1.0.1` of `electron-userland/electron-builder-binaries`. The archive holds `lib/wine/x86_64-unix/` only, without the Windows DLL directories (`x86_64-windows/`, `i386-windows/`), so Wine fails at startup: `failed to load .../lib/wine/x86_64-unix/ntdll.dll error c0000135`. Its `ntdll.so` also needs `libunwind.so.8`, which the development image does not install
- The Wine of the Windows cross image (Debian 12 `wine64` 8.0) runs only 64-bit programs, while the NSIS installer is a 32-bit program

## Impact

- Windows users keep the extracted folder themselves and start the app from it; removal is deleting the folder
- The desktop packages of Linux (AppImage) and macOS (dmg) are unaffected

## Options

- A later Wine toolset of electron-builder whose Linux archive contains the Windows DLLs, then `win.target: "nsis"` in `packages/app-desktop/electron-builder.config.ts`
- Building the Windows package on a Windows runner in CI, where NSIS needs no Wine; local builds would keep the zip

## Locations

- `packages/app-desktop/electron-builder.config.ts`: `win.target`
- `dev/desktop/package`, `dev/desktop/start-test`: the `zip` extension and the extraction of the Windows package
