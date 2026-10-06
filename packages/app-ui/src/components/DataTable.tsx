import {
  createColumnHelper,
  createSortedRowModel,
  metaHelper,
  rowSortingFeature,
  sortFns,
  tableFeatures,
  useTable,
  type ColumnDef,
  type RowData,
  type SortingState,
} from "@tanstack/react-table";
import { useCallback, useMemo } from "react";
import {
  Cell,
  Column,
  Row,
  Table,
  TableBody,
  TableHeader,
  TableLayout,
  Virtualizer,
  type ColumnProps,
  type SortDescriptor,
  type SortDirection,
} from "react-aria-components";
import ArrowDown from "~icons/lucide/arrow-down";
import ArrowUp from "~icons/lucide/arrow-up";
import ArrowUpDown from "~icons/lucide/arrow-up-down";

import { cn } from "../ui/cn";
import { sortDescriptor } from "./sortDescriptor";

interface DataColumnMeta {
  width: NonNullable<ColumnProps["defaultWidth"]>;
  minWidth?: NonNullable<ColumnProps["minWidth"]>;
  numeric?: boolean;
}

const FEATURES = tableFeatures({
  rowSortingFeature,
  sortedRowModel: createSortedRowModel(),
  sortFns,
  columnMeta: metaHelper<DataColumnMeta>(),
});

const LAYOUT = { headingHeight: 40, estimatedRowHeight: 41 };

const DEFAULT_WIDTH = "1fr";

export type DataFeatures = typeof FEATURES;

export type DataColumn<Row extends RowData> = ColumnDef<DataFeatures, Row>;

export function dataColumns<Row extends RowData>() {
  return createColumnHelper<DataFeatures, Row>();
}

export function DataTable<Data extends RowData>({
  columns,
  rows,
  rowId,
  initialSorting,
  rowClassName,
  label,
}: {
  columns: ReadonlyArray<DataColumn<Data>>;
  rows: readonly Data[];
  rowId: (row: Data) => string;
  initialSorting: SortingState;
  rowClassName?: ((row: Data) => string | undefined) | undefined;
  label: string;
}) {
  const options = useMemo(
    () => ({
      features: FEATURES,
      columns,
      data: rows,
      getRowId: rowId,
      initialState: { sorting: initialSorting },
      enableMultiSort: false,
    }),
    [columns, initialSorting, rowId, rows],
  );

  const table = useTable(options);

  const headers = table.getHeaderGroups().flatMap((group) => group.headers);
  const items = table.getRowModel().rows;
  const sort = sortDescriptor(table.state.sorting);

  const dependencies = useMemo(() => [columns, rowClassName], [columns, rowClassName]);

  const toggleSorting = useCallback(
    (descriptor: SortDescriptor) => table.getColumn(String(descriptor.column))?.toggleSorting(),
    [table],
  );

  return (
    <Virtualizer layout={TableLayout} layoutOptions={LAYOUT} shouldObserveItemSize>
      <Table
        aria-label={label}
        {...(sort === undefined ? {} : { sortDescriptor: sort })}
        onSortChange={toggleSorting}
        className="max-h-[36rem] w-full overflow-auto overscroll-contain text-xs"
      >
        <TableHeader className="bg-card border-b">
          {headers.map((header) => (
            <Column
              key={header.id}
              id={header.column.id}
              allowsSorting
              defaultWidth={header.column.columnDef.meta?.width ?? DEFAULT_WIDTH}
              minWidth={header.column.columnDef.meta?.minWidth ?? null}
              className={cn(
                "text-muted-foreground data-focus-visible:ring-ring data-hovered:text-foreground flex cursor-default items-center gap-1 px-2 font-bold outline-none data-focus-visible:ring-2 data-focus-visible:ring-inset",
                header.column.columnDef.meta?.numeric === true && "justify-end",
              )}
            >
              {({ sortDirection }) => (
                <>
                  <table.FlexRender header={header} />
                  <SortIcon direction={sortDirection} />
                </>
              )}
            </Column>
          ))}
        </TableHeader>
        <TableBody items={items} dependencies={dependencies}>
          {(row) => (
            <Row
              id={row.id}
              className={cn(
                "hover:bg-muted/50 data-focus-visible:ring-ring border-b outline-none data-focus-visible:ring-2 data-focus-visible:ring-inset",
                rowClassName?.(row.original),
              )}
            >
              {row.getAllCells().map((cell) => (
                <Cell
                  key={cell.id}
                  className={cn(
                    "data-focus-visible:ring-ring flex items-center p-2 whitespace-nowrap outline-none data-focus-visible:ring-2 data-focus-visible:ring-inset",
                    cell.column.columnDef.meta?.numeric === true && "justify-end",
                  )}
                >
                  <table.FlexRender cell={cell} />
                </Cell>
              ))}
            </Row>
          )}
        </TableBody>
      </Table>
    </Virtualizer>
  );
}

function SortIcon({ direction }: { direction: SortDirection | undefined }) {
  if (direction === "ascending") {
    return <ArrowUp aria-hidden className="size-3" />;
  }

  if (direction === "descending") {
    return <ArrowDown aria-hidden className="size-3" />;
  }

  return <ArrowUpDown aria-hidden className="size-3 opacity-40" />;
}
