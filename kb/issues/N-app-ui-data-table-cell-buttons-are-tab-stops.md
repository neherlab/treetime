# Buttons in result table cells are separate tab stops

The result tables (`DataTable` in `packages/app-ui/src/components/DataTable.tsx`) render as a React Aria grid, which lets arrow keys move between cells and into a cell's button. The cell buttons themselves (`PositionButton` in `packages/app-ui/src/runs/RecurrentTable.tsx`, `SampleButton` in `packages/app-ui/src/runs/HomoplasyTables.tsx`) are Base UI buttons with `tabIndex` 0, so the Tab key also stops at every rendered button.

## Evidence

Homoplasy run on `data/sc2/4500` in the dev app: the recurrent-mutation grid renders about 20 rows, and 19 of its buttons have a non-negative `tabIndex`. Tab walks through them one by one before it leaves the table. The WAI-ARIA grid pattern keeps one tab stop per grid and moves within it by arrow keys.

## Fix direction

Make the grid the single tab stop: let the row or cell carry the action (React Aria `onRowAction` or `onAction`), or render the cell controls with React Aria `Button`, which takes part in the grid's focus management.

> [!IMPORTANT]
> **Decision required.** A row action changes how the recurrent-mutation and sample tables look and respond to the pointer (the whole row becomes the control), while React Aria buttons keep today's look but add a second button implementation next to `packages/app-ui/src/ui/button.tsx`.
