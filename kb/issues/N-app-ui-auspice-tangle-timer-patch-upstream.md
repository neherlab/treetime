# Auspice redraws the lines between two trees from a stale timer

## Problem

When two trees are drawn side by side, as on the compare page, Auspice 3.2.0 draws lines between the tips of the same name (the tanglegram). After a resize, `componentDidUpdate()` of the tangle component removes the lines and starts a timer that draws them again after 1 s. The timer is never reset after it fires, so at the next resize the component clears the fired timer instead of removing the old lines, and the old lines stay for 1 s at their old positions. The timer is also never cleared on unmount, so a redraw that is still pending when the user leaves the page draws into an SVG that is no longer in the document.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.2.0.patch`: the timer marks itself as done when it fires, and `componentWillUnmount()` clears a pending timer.

## Code

- Redraw timer: [tree/tangle/index.js#L72-L81](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/tree/tangle/index.js#L72-L81)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.

## Proposed implementation plan

1. Open a pull request on `nextstrain/auspice` with the hunk for `src/components/tree/tangle/index.js` from `patches/auspice@3.2.0.patch`: in `Tangle.componentDidUpdate()` the timer callback sets `this.timeout = false` before `drawLines()`, and a new `Tangle.componentWillUnmount()` clears a pending `this.timeout`. Auspice 3.2.0 (tag `v3.2.0`, commit `db211cc6`), the release that the patch targets, still has the unchanged code. Check first whether upstream `master` changed it after that release.
2. Upstream test: in Auspice's Jest suite, mount `Tangle` with fake timers, change `width` twice with the timer run out between the changes, and assert that the second change removes the `.tangleLine` elements; unmount with a pending timer and assert that `drawLines` does not run. Oracle: the behavior that the code comment of `componentDidUpdate()` states (remove the lines, redraw after the resize ends).
3. After an Auspice release contains the fix, upgrade `auspice` in `package.json`, regenerate the Bun patch without the tangle hunk (procedure in [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md)), and delete this issue.
4. Validate: `just typecheck`, `just test-ts`, and a manual resize on the compare page with two trees: the old lines disappear at once at every resize, and no console error appears after leaving the page during a redraw.
