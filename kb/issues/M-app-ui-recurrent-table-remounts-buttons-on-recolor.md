# Homoplasy tables remount every button when the tree is recolored

`RecurrentTable` in `packages/app-ui/src/runs/RecurrentTable.tsx` builds its column definitions from the set of pressed mutations, and the `cell` renderer of the mutation column is a closure over that set. `AmbiguousPanel` in `packages/app-ui/src/runs/HomoplasyTables.tsx` does the same with the selected position. Every recolor of the tree therefore creates a new `cell` function.

TanStack Table renders a function `cell` as a component type: `flexRender()` in `@tanstack/react-table` calls `React.createElement(def.cell, cell.getContext())`. A new function at the same position is a different component type, so React unmounts the old cell subtree and mounts a new one. On a recolor, every `PositionButton` of the table, with its Base UI `Button` and DOM node, is destroyed and created again, although only two rows change.

## Evidence

Homoplasy run on `data/sc2/4500`, dev build in Chrome 149: after a recolor, the `<tr>` and `<td>` nodes of the recurrent-mutation table (4,470 rows) are the same nodes, and every `<button>` is a new node. With the `cell` renderer kept stable and the pressed set passed through a React context, the buttons stay mounted and the recolor takes 4.1 s instead of 4.9 s (median of three).

## Fix direction

Keep column definitions independent of selection state: define `cell` renderers as module-scope components, and let each cell or row read whether it is pressed from a React context or a per-row selector. Compute the row highlight in the row from the same source, so that a recolor re-renders only the rows whose state changes. The same rule applies to every `DataTable` column definition.
