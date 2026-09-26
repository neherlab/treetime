# The OpenAPI generator aide is pinned to a pre-release

The API server registers its routes through aide, which derives the OpenAPI document from the handler types. The workspace pins the pre-release `aide = "=0.16.0-alpha.4"` (`Cargo.toml`, used by `packages/app-server`), because the latest stable release, 0.15.1, depends on schemars 0.9, and every data type of the project derives the schemars 1.x `JsonSchema`.

## Impact

- A pre-release receives no patch releases; a fix in aide reaches the project only with a later alpha or with 0.16.0
- The API of a pre-release can change between alphas, so moving to another alpha may need code changes in `packages/app-server/src/routes.rs` and `packages/app-server/src/api/`

## Required change

Upgrade to aide 0.16.0 once it is released and at least seven days old, and check that `just generated-check` shows no change of `packages/app-contracts/openapi.json` or that each change is intended. The release is tracked in [tamasfe/aide#270](https://github.com/tamasfe/aide/issues/270).

## Evidence

- aide 0.15.1 requires `schemars = "^0.9.0"` (crates.io dependency list of the release)
- The tracking issue for the 0.16.0 release, [tamasfe/aide#270](https://github.com/tamasfe/aide/issues/270), is open
