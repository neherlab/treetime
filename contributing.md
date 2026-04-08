# Contributing to TreeTime

Thank you for your interest in contributing to TreeTime.
We welcome pull-requests that fix bugs or implement new features. 

## Bugs
If you come across a bug or unexpected behavior, please file an issue. 

## Testing
Upon pushing a commit, GitHub Actions runs tests across supported Python versions and lint checks. Tests use data from the [neherlab/treetime_examples](https://github.com/neherlab/treetime_examples) repository.

## Releasing

### Steps

1. Add a `# Unreleased` heading at the top of `changelog.md` followed by bullet points describing the changes. Commit to master (directly or via PR).
2. Run the release script from the repo root:
   ```bash
   ./release 0.x.y
   ```
3. Review the diff shown by the script, then confirm to push.
4. Watch CI complete at https://github.com/neherlab/treetime/actions. On success, the [GitHub Release](https://github.com/neherlab/treetime/releases) and the [PyPI package](https://pypi.org/project/phylo-treetime/) are created automatically.

### How it works

The [`release`](release) script:
- Validates the working tree: clean, on master, up-to-date with origin, tag absent, `# Unreleased` present in changelog.
- Bumps the version in [`treetime/__init__.py`](treetime/__init__.py).
- Renames the `# Unreleased` heading to `# 0.x.y` in [`changelog.md`](changelog.md).
- Commits (`chore: release 0.x.y`) and creates tag `v0.x.y`.
- Pushes the commit and tag after confirmation.

The tag push triggers the [release CI workflow](.github/workflows/release.yml), which:
- Runs the full test matrix and lint.
- Extracts the changelog section using [`scripts/extract-release-notes`](scripts/extract-release-notes) and creates a [GitHub Release](https://github.com/neherlab/treetime/releases).
- Builds and publishes the package to [PyPI](https://pypi.org/project/phylo-treetime/) using the `PYPI_API_TOKEN` repository secret.

   To verify the secret is present:

   ```bash
   gh secret list
   ```

## Coding conventions (loosly adhered to)

  * indentation: 4 spaces
  * docstrings: numpy style
  * variable names: snake_case
