# Package channels missing

## Summary

TreeTime v1 ships as downloads from the nightly releases and as the nightly Docker images `neherlab/treetime` on Docker Hub. Pipeline users (Nextstrain and augur, Snakemake, Nextflow) and HPC users install tools through package managers and containers, and `pip install phylo-treetime` and bioconda `treetime` deliver v0 today.

## Evidence

- Releases: `neherlab/treetime-nightly` (`dev/publish-nightly`), and the images `neherlab/treetime` of each nightly (`dev/docker/prod-images`). Apptainer on HPC pulls these images as well.
- v0 is published as `phylo-treetime` on PyPI ([packages/legacy/treetime/setup.py#L18](../../packages/legacy/treetime/setup.py#L18)) and as `treetime` on bioconda.

> [!IMPORTANT]
> **Decision required.** Which channels v1 ships through, and when. An earlier review rejected bioconda as out of scope; the v1 release reopens the question. Options, each independent:
>
> - bioconda: the usual path of bioinformatics pipelines and HPC users without root; needs the source build of [M-release-source-build-untested-outside-linux.md](M-release-source-build-untested-outside-linux.md)
> - a PyPI wheel that carries the binary, as `ruff` and `uv` do: makes `pip install` and the augur integration work; needs a plan for the `phylo-treetime` name of v0
> - Homebrew (macOS and Linux) and winget (Windows) for interactive CLI users

## Impact

- Pipelines keep v0 until a channel they use carries v1.

## Proposed implementation plan

Each channel depends on the open decision above; the options are independent, so each gets its own short plan.

- **bioconda**: after the source build works on macOS ([M-release-source-build-untested-outside-linux.md](M-release-source-build-untested-outside-linux.md)), write a `meta.yaml` recipe that builds the `treetime` binary from a release tarball with `cargo install` and the conda compilers, and submit it to `bioconda-recipes`. Name: a new package or the next major version of `treetime`, which belongs to the open decision. Validate with `bioconda-utils build` and a run of `treetime --version` in the test section
- **PyPI wheel**: build platform wheels that carry the binary with `maturin` in `bin` mode (as `ruff` does), one per target of `dev/cross/targets`, and publish them from the release workflow with trusted publishing. Settle the `phylo-treetime` name first (new major version of the same project or a new project). Validate with `pip install` of the built wheel in a clean virtual environment on each runner of `compat-*` and a run of one command
- **Homebrew and winget**: a tap formula that downloads the release asset and its checksum, updated by the release workflow; a winget manifest per release, submitted with `wingetcreate`. Validate with `brew install` and `brew test` of the formula, and `winget validate` of the manifest
- **Docs**: list each new channel on the installation page ([M-release-install-docs-and-system-requirements-missing.md](M-release-install-docs-and-system-requirements-missing.md))
