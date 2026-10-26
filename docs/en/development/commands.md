---
title: Everyday commands
description: The just recipes for routine work, and how their names are formed.
---

| Task                                             | Command                              |
| ------------------------------------------------ | ------------------------------------ |
| List the recipes                                 | `just`                               |
| Fast checks (format, clippy, TypeScript)         | `just check`                         |
| Every check; must pass before merging            | `just check-all`                     |
| Fast lint fixes and format (stage changes first) | `just fix`                           |
| Every lint fix, dylint included, and format      | `just fix-all`                       |
| Build                                            | `just build` (`just b`)              |
| Run the CLI                                      | `just run treetime ancestral --help` |
| Rust tests, optionally filtered                  | `just test-rs [filter]` (`just t`)   |
| TypeScript tests                                 | `just test-ts`                       |
| Clippy                                           | `just lint-rs` (`just l`)            |
| Clippy of the WebAssembly back end               | `just lint-wasm`, `just lint-wasm-mt` |
| Smoke test of the WebAssembly back end, in Bun   | `just test-wasm`                     |
| Smoke test of the threaded WebAssembly back end  | `just test-wasm-mt`                  |
| TypeScript lints                                 | `just lint-ts`                       |
| Format one toolchain                             | `just fmt-rs`, `fmt-ts`, `fmt-other` |
| Lints and tests of the tools in `dev/`           | `just check-dev`                     |

In the container, prefix each command with `./dev/docker/run`.

Recipe names follow one scheme:

- **Leaf recipes** run one tool in one mode. The suffix `-rs`, `-wasm`, or `-ts` names the toolchain: Rust, Rust for WebAssembly, or TypeScript, and `-wasm-mt` the threaded WebAssembly build; tools that serve one toolchain only keep their own name, such as `dylint`, `hawk`, or `knip`
- **Combined recipes** without a language suffix, such as `lint`, `test`, `fmt`, or `fix`, run the leaves of both toolchains and call no tool themselves
- **The suffix `-all`** adds the slow tools to the fast set of the same name: `fix-all` adds the dylint fixes, `lint-all` the lint libraries, hawk, and the dependency and config lints, and `check-all` the full gate

While working on one toolchain, run its leaves (`just fmt-rs`, `just l`, `just t`) and leave the slow tools to `check-all`.

The lint libraries and oxlint plugins in `dev/lints` have their own recipes, which `check` and `check-all` leave out: `just check-dev` runs `lint-dev-rs` and `test-dev-rs` for the dylint libraries, each in its own workspace with its pinned nightly, `lint-dev-ts` and `test-dev-ts` for the oxlint plugins, and `test-dev-sh` for the shell libraries in `dev/lib` (the port slots of `dev/lib/ports.sh`). Run them after a change in `dev/lints` or `dev/lib`.

`check` and `check-all` run their checks in parallel through `dev/run-checks`, keep going past failures, and list the failed checks at the end. Each check writes its output to `tmp/checks/<check>.log`. Warnings fail the checks, and a missing tool is a failure, not a skipped check. The full gate consists of groups, one CI job each; `just check-group <group>` runs one group serially, as its CI job does.

## Examples

Replace `$v` with a dataset path under `data/`, for example `flu/h3n2/20`:

```bash
just run treetime ancestral --method-anc=marginal --tree=data/$v/tree.nwk --alignment=data/$v/aln.fasta.xz --output-all=tmp/ancestral/$v
just run treetime clock --tree=data/$v/tree.nwk --metadata=data/$v/metadata.tsv --output-all=tmp/clock/$v
just run treetime timetree --tree=data/$v/tree.nwk --metadata=data/$v/metadata.tsv --alignment=data/$v/aln.fasta.xz --output-all=tmp/timetree/$v
```
