# Downloads need a choice among target triplets

## Summary

A user who downloads the CLI picks one of eight target triplets (for example `x86_64-unknown-linux-gnu` or `-musl`, `aarch64-pc-windows-gnullvm`). The static musl binaries run on every Linux distribution, so the gnu CLI downloads add a choice without a benefit for users; the gnu builds are needed only for the Node addon of the desktop app.

## Evidence

- `dev/cross/targets` lists the eight CLI targets; `dev/publish-nightly` uploads one asset per target.
- The compat jobs run the musl binaries on every distribution image of `dev/test-linux-distros-musl`, the empty `scratch` image included.

## Potential solutions

- An install script (`install.sh`, `install.ps1`) that picks the right asset for the machine.
- One Linux CLI download per architecture (musl), with the gnu binaries kept for the desktop app only.

## Proposed implementation plan

The two solutions combine: drop the Linux gnu CLI assets, then add install scripts that pick among the six remaining assets.

1. **Linux assets.** In `dev/publish-nightly`, upload the CLI only for the musl Linux targets; keep building the gnu targets of `dev/cross/targets` for the `app-napi` addon of the desktop app. Update the download list in the docs (`docs/en/`) to one Linux file per architecture
2. **Install scripts.** Add `install.sh` (Linux, macOS: map `uname -s` and `uname -m` to the musl or darwin triplet) and `install.ps1` (Windows: map `$env:PROCESSOR_ARCHITECTURE` to the gnullvm triplet), upload them as release assets, download the matching asset, check its checksum, and install into a user directory
3. **Tests.** Factor the platform-to-triplet mapping of `install.sh` into a function and test it with `bats` or a shell test recipe over the table of `uname` outputs (oracle: `dev/cross/targets`); run `install.sh` in the musl distribution images of `dev/test-linux-distros-musl` against a nightly release
4. **Validate.** `just lint-shell`, `just fmt-other`, `just lint-workflows` if a workflow changes
