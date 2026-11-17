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
import { cn } from "cn";
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
  type Key,
  type Selection,
  type SortDescriptor,
  type SortDirection,
} from "react-aria-components";
import ArrowDown from "~icons/lucide/arrow-down";
import ArrowUp from "~icons/lucide/arrow-up";
import ArrowUpDown from "~icons/lucide/arrow-up-down";

import { sortDescriptor } from "./sort-descriptor";

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
  selection,
}: {
  columns: ReadonlyArray<DataColumn<Data>>;
  rows: readonly Data[];
  rowId: (row: Data) => string;
  initialSorting: SortingState;
  rowClassName?: ((row: Data) => string | undefined) | undefined;
  label: string;
  selection?: RowSelection | undefined;
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

  const dependencies = useMemo(() => [columns, rowClassName], [columns, rowClassName]);

  const selectionProps = useMemo(() => tableSelection(selection), [selection]);

  const toggleSorting = useCallback(
    (descriptor: SortDescriptor) => table.getColumn(String(descriptor.column))?.toggleSorting(),
    [table],
  );

  return useMemo(() => {
    const headers = table.getHeaderGroups().flatMap((group) => group.headers);
    const sort = sortDescriptor(table.state.sorting);

    return (
      <Virtualizer layout={TableLayout} layoutOptions={LAYOUT} shouldObserveItemSize>
        <Table
          aria-label={label}
          {...(sort === undefined ? {} : { sortDescriptor: sort })}
          onSortChange={toggleSorting}
          {...selectionProps}
          className="max-h-[36rem] w-full overflow-auto overscroll-contain text-xs"
        >
          <TableHeader className="bg-card size-full border-b">
            {headers.map((header, index) => (
              <Column
                key={header.id}
                id={header.column.id}
                isRowHeader={index === 0}
                allowsSorting
                defaultWidth={header.column.columnDef.meta?.width ?? DEFAULT_WIDTH}
                minWidth={header.column.columnDef.meta?.minWidth ?? null}
                className={cn(
                  "text-muted-foreground data-focus-visible:ring-ring data-hovered:text-foreground flex h-full items-center gap-1 px-2 font-bold outline-none data-focus-visible:ring-2 data-focus-visible:ring-inset",
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
          <TableBody items={table.getRowModel().rows} dependencies={dependencies}>
            {(row) => (
              <Row
                id={row.id}
                className={cn(
                  "hover:not-data-selected:bg-muted/50 data-selected:bg-selected data-focus-visible:ring-ring size-full border-b outline-none data-focus-visible:ring-2 data-focus-visible:ring-inset",
                  rowClassName?.(row.original),
                )}
              >
                {row.getAllCells().map((cell) => (
                  <Cell
                    key={cell.id}
                    className={cn(
                      "data-focus-visible:ring-ring flex h-full items-center overflow-hidden p-2 whitespace-nowrap outline-none data-focus-visible:ring-2 data-focus-visible:ring-inset",
                      cell.column.columnDef.meta?.numeric === true && "justify-end",
                    )}
                  >
                    <div className="min-w-0 truncate">
                      <table.FlexRender cell={cell} />
                    </div>
                  </Cell>
                ))}
              </Row>
            )}
          </TableBody>
        </Table>
      </Virtualizer>
    );
  }, [dependencies, label, rowClassName, selectionProps, table, toggleSorting]);
}

export interface RowSelection {
  selected: string | undefined;
  onChange: (id: string) => void;
}

export function CellText({ text }: { text: string }) {
  return <span title={text}>{text}</span>;
}

function tableSelection(selection: RowSelection | undefined) {
  if (selection === undefined) {
    return {};
  }

  const { selected, onChange } = selection;

  return {
    selectionMode: "single",
    selectionBehavior: "replace",
    disallowEmptySelection: true,
    selectedKeys: selected === undefined ? [] : [selected],
    onSelectionChange: (keys: Selection) => {
      const key = selectedKey(keys);

      if (key !== undefined) {
        onChange(String(key));
      }
    },
  } as const;
}

function selectedKey(keys: Selection): Key | undefined {
  return keys === "all" ? undefined : keys.values().next().value;
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
