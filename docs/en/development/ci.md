---
title: Continuous integration
description: The CI workflows and their container images.
---

The workflows in `.github/workflows/` run the same `just` recipes and scripts as a developer, in the same containers. Deployments never wait for checks or tests: a check reports a problem, and a failed target, image, or test never holds back the others.

- **Pull request**: `cli.yml` runs the check groups of `just check-all` in parallel jobs, except `hawk`, which needs a second full compilation; run it locally. `web.yml` and `docs.yml` build the web2 app and the docs site and publish them to the preview site of the pull request (see [Web app deploys](/development/releases#web-app-deploys) and [Docs deploys](/development/releases#docs-deploys)). Pull requests from forks build without an upload, because they get no secrets
- **Push to any branch**: `web.yml` and `docs.yml` publish the preview sites of the commit, of the branch, and of each open pull request of the branch. A pull request starts a run of its own only when it is opened or reopened
- **Push to `rust`**: the check groups, the production web app and docs site, and the nightly: `nightly.yml` builds the CLI for every target once, releases it, publishes the Docker images, and tests the binaries (see [Nightly releases](/development/releases#nightly-releases))
- **Manual run of `CI` on a branch** (`gh workflow run cli.yml --ref <branch>`; the Actions page shows no button, because the workflow is not on `master`): the check groups, and the CLI build of every target in dev mode with the compatibility tests of `cli-compat.yml`, to test a change of the build before it lands
- **Daily at 04:00 UTC**: `schedule-nightly.yml` on `master` calls the nightly, which skips when the latest nightly already has the tip of `rust`

`master` holds TreeTime 0.x and its own workflows. No workflow of `rust` runs on `master` or changes it; `schedule-nightly.yml` is there only because GitHub runs scheduled workflows from the default branch.

The CI jobs pull the container images from Docker Hub by the hash of their build inputs, and build them when the inputs changed. Pull requests never publish images; pushes to `rust`, the nightly, and manual runs publish them, with zstd-compressed layers, which compress and unpack faster than gzip and need Docker 23 or later to pull. The layers of each image go from the inputs that change least to those that change most, so a tool update in `.config/mise.toml` rebuilds only the last layers.
