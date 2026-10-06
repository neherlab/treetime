# Result tables render every row

`DataTable` in `packages/app-ui/src/components/DataTable.tsx` renders one table row per data row and scrolls them inside a fixed-height box. Large runs produce tables with thousands of rows, and every re-render of such a table renders all of them.

## Evidence

Homoplasy run on `data/sc2/4500`: the recurrent-mutation table has 4,470 rows and the sample table 4,717. Coloring the tree by another mutation re-renders the recurrent-mutation table, whose pressed button and row highlight change. The ancestral view has the same structure: its branch table lists every branch with mutations, 9,000 rows on `sc2/4500`.

Dev build in Chrome 149, median of three recolors after a fresh load:

| Variant                                                                                                                                                | Recolor |
| ------------------------------------------------------------------------------------------------------------------------------------------------------ | ------- |
| Current code                                                                                                                                           | 4.9 s   |
| Mutation cells kept mounted ([M-app-ui-recurrent-table-remounts-buttons-on-recolor.md](M-app-ui-recurrent-table-remounts-buttons-on-recolor.md) fixed) | 4.1 s   |
| Recurrent-mutation table limited to 40 rows                                                                                                            | 2.7 s   |

The full table therefore costs about 2.2 s of the update in the dev build. The remaining 2.7 s is outside the table (Auspice and the other panels) and does not change with the number of rows.

Dev timings overstate production cost. The CPU profile of a recolor is dominated by development-only work: owner stacks through `console.createTask`, React performance tracks through `performance.measure`, the React DevTools hook, and the second render of every component under `StrictMode` (`packages/app-web/src/main.tsx`). Without a DOM, the render phase of the 4,470-row recurrent-mutation table (`react-dom/server` `renderToString`, real rows) takes 210 ms with production React and 556 ms with development React; 40 rows take 1.7 ms. This leaves out DOM creation and effects, so it is a lower bound of one full render in production.

## Fix direction

Render only the rows inside the scroll box, keeping the sticky header, sorting, and the column layout of `DataTable`.

- **Prerequisites**: stable table inputs ([N-app-ui-data-table-copies-rows-and-columns.md](N-app-ui-data-table-copies-rows-and-columns.md)), and one scroll element, which the sticky header also needs ([M-app-ui-data-table-sticky-header-scrolls-away.md](M-app-ui-data-table-sticky-header-scrolls-away.md))
- **Required behavior**: measured heights for rows that wrap (the mutation lists of the sample table); the row that holds focus stays mounted when it scrolls out of the rendered range, or the focus is lost; `aria-rowcount` on the table and a 1-based `aria-rowindex` on every rendered row
- **Accepted costs**: browser find-in-page, printing, and screen readers reach only the rendered rows

> [!IMPORTANT]
> **Decision required.** The library that renders the virtualized table is open:
>
> - **React Aria Components `Table` with `Virtualizer` and `TableLayout`**: virtualization, sticky header, measured row heights (`estimatedRowHeight`), the focused row kept mounted, ARIA row counts, one tab stop with arrow-key navigation, and row actions are part of the library. Rows are `div` elements with `role="grid"`, and column widths come from the layout (`defaultWidth`, default `1fr`) instead of cell content. The project rules name React Aria Components the primary UI library, but no manifest lists it yet; it replaces the rendering of `DataTable` and its call sites
> - **react-virtuoso `TableVirtuoso`**: native `<table>` with spacer rows, a sticky `<thead>` whose height it measures, and automatic row measurement; an official example combines it with TanStack Table, so the row model and header of `DataTable` stay. Focus persistence and ARIA row counts are not part of the library. It adds a second library for virtual lists next to `@tanstack/react-virtual` (`packages/app-ui/src/runs/LogLines.tsx`)
> - **`@tanstack/react-virtual`**, the installed owner of virtual lists: computes the rendered range and measures rows; spacer rows or a grid layout, the header offset (`scrollMargin`), focus persistence (`rangeExtractor`), and ARIA row counts are project code
