# Auspice tree keeps drawing after it unmounts

## Problem

Auspice 3.2.0's tree component removes its drawing when it unmounts, but the tree objects keep their pending drawing steps: staged transitions of an update, and the timers that draw tip labels. When a result view closes during an update (for example, the user opens another run right after a change of the coloring), these steps run against a drawing that is no longer in the document.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.2.0.patch`, registered in `patchedDependencies` of the root `package.json`: `componentWillUnmount()` calls a `dispose()` step on each tree object before the drawing is removed, which stops the staged transitions and the tip label timers.

## Code

- `TreeComponent.componentWillUnmount()`: removes the drawing ([tree/tree.tsx#L81-L86](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/tree/tree.tsx#L81-L86))

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice). After upgrading to an Auspice release that contains it, remove the change from the patch. The patch file also holds the changes in the other `N-app-ui-auspice-*-patch-upstream.md` issues. When the patch has no changes left, delete the patch file and its `patchedDependencies` entry. Bun refuses to install when the patch no longer applies to a new Auspice version, so an upgrade without this step fails visibly.

Keep each fix of an Auspice defect free of TreeTime specifics when editing it, so that an upstream contribution can take its code unchanged and remove it from the patch.

## Editing the patch

The patch is a Git diff against the published Auspice package. Never write it with `bun patch --commit`, because Bun adds its `.bun-tag-*` marker files to the patch. Edit it against a copy of the published package instead:

1. Run `npm pack auspice@<version>` in a scratch directory, and extract the archive twice: once as the base, once for editing
2. In the copy for editing, apply the current patch with `patch -p1 < <treetime>/patches/auspice@<version>.patch`, and edit the files
3. Run `git init`, `git add -A` and `git commit` in the base copy, copy the edited `src/` over it, and write the patch with `git diff --full-index > auspice@<version>.patch`
4. Replace the patch file, delete `node_modules/auspice`, `node_modules/.bun/auspice@*` and the patched copies of Auspice in the Bun install cache (`auspice@<version>@@@1_patch_hash=*` in `.cache/bun` of the checkout), and run `./dev/docker/run bun install`. Bun can otherwise install a patched copy from an earlier version of the patch file. Check that the files in `packages/app-ui/node_modules/auspice` contain the change
5. Delete Vite's dependency cache (`node_modules/.vite` of the app packages) and restart the dev server, because Vite bundles Auspice again only when `bun.lock` or the file list of `patches/` changes, and an edit inside the patch file changes neither

## Proposed implementation plan

1. Open a pull request on `nextstrain/auspice` with the `dispose()` hunks of `patches/auspice@3.2.0.patch`: `renderers.dispose()` sets `disposed`, stops the staged transitions and clears the tip label timers, `PhyloTree.prototype.dispose` exposes it, the drawing steps return early on a disposed tree, and `TreeComponent.componentWillUnmount()` calls `dispose()` on each tree object before it removes the drawing. Auspice 3.2.0 (tag `v3.2.0`, commit `db211cc6`), the release that the patch targets, has no `dispose` step. Check first whether upstream `master` added one after that release.
2. Upstream test: with fake timers, start an update of a `PhyloTree` that stages transitions and tip label timers, call `dispose()`, run all timers, and assert that no element is added to or changed in the SVG. Oracle: a disposed tree must not draw (the contract of `dispose()`).
3. After an Auspice release contains the fix, upgrade `auspice` in `package.json`, regenerate the patch without these hunks with the steps of "Editing the patch", and delete this issue. Move the "Editing the patch" section to the remaining `N-app-ui-auspice-*-patch-upstream.md` issue or to the developer guide in `docs/en/development/` first, because the other issues link to it here.
4. Validate: `just typecheck`, `just test-ts`, and in the browser: change the coloring of a large tree and open another run at once; the console shows no error from the closed view.
