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
