---
title: Dependencies
description: The release age rule and how dependencies are pinned and updated.
---

Every dependency release must be at least seven days old before the project adopts it:

- Bun refuses younger npm packages (`minimumReleaseAge` in `bunfig.toml`)
- Cargo has the age check as an unstable feature until Rust 1.100. `just deps-update` and `just deps-upgrade` enable it for their resolution, and `just deps-age` fails on a lockfile entry younger than seven days
- Python packages: `exclude-newer` in `pixi.toml` makes `pixi lock` skip younger releases
- Tools in `.config/mise.toml`: `minimum_release_age` makes `just tools-outdated` list only newer releases that are at least seven days old. mise does not filter an exact version, so check the date of a version you type by hand
- Base images in `dev/docker/` follow the same rule by hand
- Releases of the [binarylandia](https://github.com/binarylandia) GitHub organization (the prebuilt OpenBLAS, liblzma, and libzstd archives, the compilers, the patched `wasm-bindgen` CLI) are the project's own builds and need no waiting time; `minimum_release_age_excludes` in `.config/mise.toml` lists their repositories. The sources that such a release builds must still be at least seven days old
- mise itself: the images and CI install the version in `dev/docker/files/mise-version`, verified by its line in `dev/docker/files/checksums`. `min_version` in `.config/mise.toml` is the oldest mise that reads the configuration; raise it when the configuration needs a newer mise
- After changing a tool in `.config/mise.toml`, `just tools-lock <tool>` locks only that tool, and after changing an entry of the images in `.config/mise.image.toml`, `just tools-lock-image <tool>` locks it in `.config/mise.image.lock`. Locking calls the GitHub API, which allows 60 anonymous requests per hour. With `MISE_GITHUB_TOKEN` set, mise authenticates instead, on the host and in the container alike (`dev/docker/run` forwards it), for example `MISE_GITHUB_TOKEN="$(gh auth token)" ./dev/docker/run just tools-lock <tool>`

Rust dependencies are pinned exactly in the workspace `Cargo.toml`, JavaScript dependencies exactly in the manifests, with shared packages in the Bun catalog of the root `package.json`. The React packages stay on 19.x, because Auspice runs in-process and each React major needs its own check of the Auspice modules TreeTime imports; `just lint-ts` enforces this.

`just deps-update` changes `Cargo.lock` in any checkout: after a new dependency in `Cargo.toml`, `just deps-update --workspace` adds it to the lockfile without updating the other entries. The upgrade recipes `just deps-upgrade` and `just deps-upgrade-ts` run in the main checkout only. `just audit` checks both dependency graphs against the security advisory databases.
