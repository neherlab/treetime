---
title: Lint settings
description: Where the project keeps the settings of its lints, and how to suppress a lint.
---

The custom dylint library and the oxlint plugin of `dev/lints/oxlint/` are both named `custom`, so a suppression reads `#[cfg_attr(dylint_lib = "custom", expect(no_comments, reason = ".."))]` or `// oxlint-disable-next-line custom/no-vague-identifiers -- reason`.

Clippy rejects `#[allow]` on an item: suppress a lint on the narrowest item with `#[expect(<lint>, reason = "..")]`, which fails once the lint no longer fires, so a stale suppression cannot stay.

The project settings of the lints are in these files:

- `justfile` variables: `dylint_rustflags` (the lint levels of the custom library), the check groups `checks_*` with `check_fast` and `check_full`, `lint_fast`, `lint_full`, `public_api_crates`, `dockerfiles`, and `react_pin_reason`
- `dylint.toml`: the settings of the dylint libraries, such as the render sources of `no_comments`, the entry points of `forbidden-reach`, the error style of `proper_error_type`, and the helper macros of `prefer_error_macros`
- `oxlint.config.ts`: a call of `projectConfig` from `dev/lints/oxlint/config.ts` with the package layout, the web scopes, the import and property restrictions, and the contracts package
- `clippy.toml` (thresholds and banned types and functions; `dev/lints/clippy.toml` replaces it for the lint libraries, with the same thresholds and without the bans), `oxfmt.config.ts`, `taplo.toml`, `.config/knip.json`, `.config/deny.toml`, `.config/hawk.toml`, `.config/jscpd.json`, and `.config/mordant-baseline.toml`

Rust comments are banned in TreeTime: every lint of the custom library is on, so `no_comments` keeps only the doc comments that a render source of `dylint.toml` carries into user-facing output, and comments that start with a marker such as `SAFETY:` or `TODO`. TypeScript comments are banned by the oxlint base configuration `dev/lints/oxlint/config.ts`, apart from tool directives.
