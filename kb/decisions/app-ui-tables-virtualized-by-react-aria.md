# Result tables are virtualized by React Aria

The shared result table `DataTable` (`packages/app-ui/src/components/DataTable.tsx`) renders through React Aria Components `Table` inside a `Virtualizer` with `TableLayout`. Only the rows in or near the visible part of the scroll box are in the page. TanStack Table stays the row model: it sorts the rows and provides the column definitions and cell renderers.

## Behavior

- **Rendered rows**: a table holds the header and about 20 rows of the visible area; the ARIA row count on the table and the row index on each row give the full size to assistive technology
- **Row heights**: a single-line row is 41 px; React Aria measures every rendered cell and sizes the row to its tallest cell, so wrapped mutation lists get their full height. It measures again when a row's height changes, for example when the panel width changes
- **Columns**: each column declares its width and minimum width in its TanStack column `meta` (`DataColumnMeta`), because the virtualized layout sizes columns by declaration instead of content. Text that exceeds its column is cut with an ellipsis; name columns show the full name as a tooltip through `CellText`
- **Sorting**: header clicks run TanStack's sort toggle. The first click sorts numeric columns descending and text columns ascending, and the third click clears the sort. React Aria only shows the state through its sort descriptor, because its own toggle always starts ascending and never clears
- **Keyboard**: arrow keys move focus between cells and into the buttons of a cell; React Aria keeps the row that holds focus rendered when it scrolls out of view
- **Accepted costs**: browser find-in-page, printing, and screen readers reach only the rendered rows

## Reason

Results of large runs have tables with thousands of rows: 4,470 recurrent mutations and 4,717 samples on the homoplasy page of `data/sc2/4500`, 6,537 branches with mutations on its ancestral page. Rendering every row made each update of such a table cost about 2 s in the dev build. With virtualization, a recolor of the tree on that homoplasy page dropped from 5.45 s to 2.61 s in the dev build; the rest of that time is spent outside the tables.

`@tanstack/react-virtual`, the library for virtual lists (`packages/app-ui/src/runs/LogLines.tsx`), only computes the visible range and measures items. A table on it needs project code for the positioning of rows in table layout, the sticky header offset, keeping the focused row rendered, and the ARIA row counts. React Aria provides all of these for tables, so result tables use it as a second virtualization library, and log lines stay on TanStack Virtual.

## Implementation

- `packages/app-ui/src/components/DataTable.tsx`: `DataTable`, `DataColumnMeta`, `CellText`
- `packages/app-ui/src/components/sortDescriptor.ts`: the React Aria sort descriptor of a TanStack sorting state
- `packages/app-ui/src/runs/RecurrentTable.tsx`, `packages/app-ui/src/runs/HomoplasyTables.tsx`: the pressed mutation and the selected position reach the cells through a React context, so a recolor re-renders the visible cells without remounting them; the row of a pressed button is highlighted through CSS (`has-aria-pressed:bg-accent`)
