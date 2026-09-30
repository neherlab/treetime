# taplo installs without checksum verification

## Summary

`.config/mise.lock` records a URL but no sha256 for taplo on all four platforms. mise verifies a downloaded archive only against a recorded checksum, so the host (`mise install`) and the development image (`dev/docker/native.dockerfile`) install taplo without checking it. Every other tool in the lockfile has a checksum.

## Cause

The taplo releases publish no checksum files, and the standard `mise lock` path (`just tools-lock`) records a checksum only when the backend supplies one. mise 2026.9.10 computes a missing checksum from the downloaded artifact only in the generate lockfile mode (`lockfile_mode = "generate"`, `complete_artifact_checksums()` in mise `src/lockfile/generate.rs`). protoc had the same gap and got its checksums this way.

## Evidence

- `.config/mise.lock`: the four `[tools.taplo."platforms.*"]` tables have `url` and `url_api` but no `checksum`
- A protoc lock entry with a wrong checksum fails `mise install` with `Expected: sha256:... Actual: sha256:...`, so a recorded checksum is enforced

## Resolution options

- Run `MISE_LOCKFILE_MODE=generate ./dev/docker/run just tools-lock taplo` and commit the added checksums
- Set `MISE_LOCKFILE_MODE=generate` in the `tools-lock` recipe, so every future lock computes missing checksums. Check first that generate mode leaves the entries of other tools unchanged
