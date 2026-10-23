# Auspice tree keeps drawing after it unmounts

## Problem

Auspice 3.1.0's tree component removes its drawing when it unmounts, but the tree objects keep their pending drawing steps: staged transitions of an update, and the timers that draw tip labels. When a result view closes during an update (for example, the user opens another run right after a change of the coloring), these steps run against a drawing that is no longer in the document.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.1.0.patch`, registered in `patchedDependencies` of the root `package.json`: `componentWillUnmount()` calls a `dispose()` step on each tree object before the drawing is removed, which stops the staged transitions and the tip label timers.

## Code

- `TreeComponent.componentWillUnmount()`: removes the drawing ([tree/tree.tsx#L81-L86](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/tree/tree.tsx#L81-L86))

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice). After upgrading to an Auspice release that contains it, remove the change from the patch. The patch file also holds the changes in the other `N-app-ui-auspice-*-patch-upstream.md` issues. When the patch has no changes left, delete the patch file and its `patchedDependencies` entry. Bun refuses to install when the patch no longer applies to a new Auspice version, so an upgrade without this step fails visibly.

The patch carries the same code as the Auspice patches of Nextclade, treeknit-rs and auspice-vscode for every fix that more than one of these viewers needs. Keep these parts identical when editing them, so that one upstream contribution removes the change from all four patches.

## Editing the patch

The patch is a Git diff against the published Auspice package. Never write it with `bun patch --commit`, because Bun adds its `.bun-tag-*` marker files to the patch. Edit it against a copy of the published package instead:

1. Run `npm pack auspice@<version>` in a scratch directory, and extract the archive twice: once as the base, once for editing
2. In the copy for editing, apply the current patch with `patch -p1 < <treetime>/patches/auspice@<version>.patch`, and edit the files
3. Run `git init`, `git add -A` and `git commit` in the base copy, copy the edited `src/` over it, and write the patch with `git diff --full-index > auspice@<version>.patch`
4. Replace the patch file, delete `node_modules/auspice`, `node_modules/.bun/auspice@*` and the patched copies of Auspice in the Bun install cache (`auspice@<version>@@@1_patch_hash=*` in `.cache/bun` of the checkout), and run `./dev/docker/run bun install`. Bun can otherwise install a patched copy from an earlier version of the patch file. Check that the files in `packages/app-ui/node_modules/auspice` contain the change
5. Delete Vite's dependency cache (`node_modules/.vite` of the app packages) and restart the dev server, because Vite bundles Auspice again only when `bun.lock` or the file list of `patches/` changes, and an edit inside the patch file changes neither
