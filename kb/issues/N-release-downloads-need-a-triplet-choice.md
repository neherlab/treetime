# Downloads need a choice among target triplets

## Summary

A user who downloads the CLI picks one of eight target triplets (for example `x86_64-unknown-linux-gnu` or `-musl`, `aarch64-pc-windows-gnullvm`). The static musl binaries run on every Linux distribution, so the gnu CLI downloads add a choice without a benefit for users; the gnu builds are needed only for the Node addon of the desktop app.

## Evidence

- `dev/cross/targets` lists the eight CLI targets; `dev/publish-nightly` uploads one asset per target.
- The compat jobs run the musl binaries on every distribution image of `dev/test-linux-distros-musl`, the empty `scratch` image included.

## Potential solutions

- An install script (`install.sh`, `install.ps1`) that picks the right asset for the machine.
- One Linux CLI download per architecture (musl), with the gnu binaries kept for the desktop app only.
