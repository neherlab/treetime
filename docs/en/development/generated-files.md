---
title: Generated files
description: The committed files that `just gen` writes, and their checks.
---

The JSON schemas, the OpenAPI document, its TypeScript client, the TypeScript declarations of the Node addon and of the WebAssembly back end, the CLI reference page of the docs site, the compiled UI messages, and the index of the example data `data/index.json` are generated and committed. `just gen` regenerates them, and `just generated-check`, part of `check-all`, fails when a committed copy is stale. Never edit them by hand.

`just data-reindex` (group `data-index` of `dev/generated`) writes `data/index.json` from the files that `git ls-files data` lists, through the `generate-data-index` binary of `packages/app-cli`; run it after adding, removing, or changing a file in `data/`, because `just generated-check` fails on a stale index and the apps offer only the files that the index lists.

`just fixtures-check` checks the reference and golden-master test fixtures against `dev/registry/reference-files.toml`.
