# Auspice redraws the lines between two trees from a stale timer

## Problem

When two trees are drawn side by side, as on the compare page, Auspice 3.1.0 draws lines between the tips of the same name (the tanglegram). After a resize, `componentDidUpdate()` of the tangle component removes the lines and starts a timer that draws them again after 1 s. The timer is never reset after it fires, so at the next resize the component clears the fired timer instead of removing the old lines, and the old lines stay for 1 s at their old positions. The timer is also never cleared on unmount, so a redraw that is still pending when the user leaves the page draws into an SVG that is no longer in the document.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.1.0.patch`: the timer marks itself as done when it fires, and `componentWillUnmount()` clears a pending timer.

## Code

- Redraw timer: [tree/tangle/index.js#L72-L81](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/tree/tangle/index.js#L72-L81)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
