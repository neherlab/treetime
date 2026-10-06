# Result tables render every row

`DataTable` in `packages/app-ui/src/components/DataTable.tsx` renders one table row per data row and scrolls them inside a fixed-height box. Large runs produce tables with thousands of rows, and every re-render of such a table rebuilds all of them.

## Evidence

Dev build in Chrome, homoplasy run on `data/sc2/4500`: the recurrent-mutation table has 4,470 rows and the sample table 4,717. Coloring the tree by another mutation re-renders the recurrent-mutation table, whose pressed button and row highlight change. The update is one main-thread task of 3.7 s; with the table limited to 50 rows the same update takes 1.7 s, so the table accounts for about 2 s. The ancestral view has the same structure: its branch table lists every branch with mutations, 9,000 rows on `sc2/4500`.

## Fix direction

Render only the visible rows with `@tanstack/react-virtual`, the owner of virtual lists in the tech stack, keeping the sticky header, sorting, and the table semantics of `DataTable`. Rows that wrap (the mutation lists of the sample table) need measured row heights.

> [!IMPORTANT]
> **Investigation required.** Measure the same update on the production build (`just serve`), where React runs without development checks, to set the priority of this change.
