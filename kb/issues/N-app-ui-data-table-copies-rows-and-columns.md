# `DataTable` rebuilds its row model on every render

`DataTable` in `packages/app-ui/src/components/DataTable.tsx` passes `columns: [...columns]` and `data: [...rows]` to `useTable()`. Both are new arrays on every render.

TanStack Table 9 memoizes its models by reference: the core row model depends on `options.data`, the sorted row model on the core row model, and the cells of a row on the column objects. New arrays on every render therefore rebuild every `Row` and every cell and sort all rows again, whenever the table re-renders for any reason. The TanStack data guide states that a new `data` reference rebuilds "every row and cell object" on every render. The copies also defeat the fine-grained subscriptions of TanStack Table 9 (`table.Subscribe` with a selector), which assume stable inputs.

The copies exist because the props are read-only arrays and the `useTable` options are mutable array types.

## Fix direction

Pass the caller's arrays to `useTable()` without copying, with prop types that `useTable()` accepts, so the row model is rebuilt only when the rows or the columns change.

> [!IMPORTANT]
> **Investigation required.** The share of this rebuild in the 2.2 s that the 4,470-row recurrent-mutation table adds to a recolor ([M-app-ui-data-table-renders-every-row.md](M-app-ui-data-table-renders-every-row.md)) is not measured separately.
