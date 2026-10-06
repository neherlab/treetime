# The header of result tables scrolls out of view

`DataTable` in `packages/app-ui/src/components/DataTable.tsx` scrolls its rows in a `max-h-[36rem] overflow-auto` box and gives the `<thead>` `position: sticky; top: 0`, so the column headers should stay visible. They scroll away with the rows.

The cause is the `Table` component in `packages/app-ui/src/ui/table.tsx`: it wraps the `<table>` in a `div` with `overflow-x-auto`. When one overflow axis is not `visible`, CSS computes the other axis as `auto` too, so this `div` is a scroll container. A sticky element sticks to its nearest scroll container, which is this `div`. The `div` is as tall as the table and never scrolls vertically, so the header never sticks.

## Evidence

Homoplasy run on `data/sc2/4500`, recurrent-mutation table in Chrome 149: the computed overflow of the wrapper is `auto` on both axes. After the outer box scrolls by 2,000 px, the top of the `<thead>` is 2,000 px above the top of the box.

## Fix direction

Make the outer box of `DataTable` the only scroll container of the table, for example by rendering the `<table>` without the `overflow-x-auto` wrapper in `DataTable`; the outer box already scrolls on both axes. Row virtualization ([M-app-ui-data-table-renders-every-row.md](M-app-ui-data-table-renders-every-row.md)) needs the same single scroll element.
