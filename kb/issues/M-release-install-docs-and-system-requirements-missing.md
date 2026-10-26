# Install documentation and system requirements missing

## Summary

The `rust` branch has no root README and no user-facing page that says how to install TreeTime v1 or which systems it runs on. Users learn the CPU and operating system floors only when a binary or the desktop app refuses to start, and nothing points first-time users to the web app, which needs no installation.

## Evidence

- The repository root holds no `README.md`; the docs site holds a command overview with placeholder text (`docs/en/guide/commands.md`) and the generated command-line reference (`packages/app-docs/src/generated/en/reference/cli.md`); none of them covers installation.
- The floors exist only in the developer guide: [docs/en/development/releases.md#L22](../../docs/en/development/releases.md#L22) lists the CPU baselines of the release binaries (x86-64-v2, Armv8.0 on Linux aarch64, Armv8.1 LSE on Windows aarch64, Apple M1).
- The desktop app adds its own floors through Electron 44.4.1 ([packages/app-desktop/package.json#L25](../../packages/app-desktop/package.json#L25)), which are higher than those of the CLI for macOS and Linux.

> [!IMPORTANT]
> **Investigation required.** The macOS, Windows, and Linux (glibc) minimums of Electron 44 are not confirmed. Read them from the Electron 44 release notes and the Chromium platform support pages before the page states them.

## Impact

- A user on an unsupported system downloads a binary and gets an error message or a failed start instead of a statement before the download.
- Users who only want to try TreeTime do not learn that the web app runs in the browser without installation.

## Potential solutions

- A root `README.md` and an installation page with, for each channel (CLI download, desktop app, web app, later package managers), the steps and the system requirements: CPU baselines, oldest operating systems (glibc 2.17, Windows 10, macOS 10.12 x86_64 and 11.0 arm64), and the Electron floors of the desktop app.
- Point first-time users to the web app and name its limits (4 GiB of WebAssembly memory, one browser tab).
