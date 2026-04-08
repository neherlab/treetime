# Contributing to TreeTime

Thank you for your interest in contributing to TreeTime.
We welcome pull-requests that fix bugs or implement new features. 

## Bugs
If you come across a bug or unexpected behavior, please file an issue. 

## Testing
Upon pushing a commit, GitHub Actions runs tests across supported Python versions and lint checks. Tests use data from the [neherlab/treetime_examples](https://github.com/neherlab/treetime_examples) repository.

## Releasing

Releases are handled by the `release` script at the repo root. Before running it, add a `# Unreleased: <title>` section at the top of `changelog.md` with bullet points describing the changes.

```
./release 0.x.y
```

The script will:
1. Validate the working tree, branch, and changelog
2. Set the version in `treetime/__init__.py` and rename `# Unreleased` to `# 0.x.y` in `changelog.md`
3. Show a diff and commit `chore: release 0.x.y`
4. Create tag `v0.x.y`
5. Ask for confirmation, then push the commit and tag

CI then runs on the tag, and on success:
- Creates a GitHub Release with the changelog section as release notes
- Publishes to PyPI

Requires a `PYPI_API_TOKEN` secret in the repository. To check it's there:

```
gh secret list
```

## Coding conventions (loosly adhered to)

  * indentation: 4 spaces
  * docstrings: numpy style
  * variable names: snake_case
