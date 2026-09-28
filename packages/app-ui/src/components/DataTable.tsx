import {
  createColumnHelper,
  createSortedRowModel,
  rowSortingFeature,
  sortFns,
  tableFeatures,
  useTable,
  type ColumnDef,
  type RowData,
  type Header,
  type SortingState,
} from "@tanstack/react-table";
import { ArrowDown, ArrowUp, ArrowUpDown } from "lucide-react";

import { Button } from "../ui/button";
import { cn } from "../ui/cn";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "../ui/table";

const FEATURES = tableFeatures({ rowSortingFeature, sortedRowModel: createSortedRowModel(), sortFns });

export type DataFeatures = typeof FEATURES;

export type DataColumn<Row extends RowData> = ColumnDef<DataFeatures, Row>;

export function dataColumns<Row extends RowData>() {
  return createColumnHelper<DataFeatures, Row>();
}

export function DataTable<Row extends RowData>({
  columns,
  rows,
  rowId,
  initialSorting,
  numeric,
  rowClassName,
  label,
}: {
  columns: ReadonlyArray<DataColumn<Row>>;
  rows: readonly Row[];
  rowId: (row: Row) => string;
  initialSorting: SortingState;
  numeric: ReadonlySet<string>;
  rowClassName?: ((row: Row) => string | undefined) | undefined;
  label: string;
}) {
  const table = useTable({
    features: FEATURES,
    columns: [...columns],
    data: [...rows],
    getRowId: rowId,
    initialState: { sorting: initialSorting },
    enableMultiSort: false,
  });

  return (
    <div className="max-h-[36rem] overflow-auto overscroll-contain">
      <Table aria-label={label} className="text-xs">
        <TableHeader className="bg-card sticky top-0 z-10">
          {table.getHeaderGroups().map((group) => (
            <TableRow key={group.id}>
              {group.headers.map((header) => (
                <TableHead
                  key={header.id}
                  aria-sort={ariaSort(header)}
                  className={cn(numeric.has(header.column.id) && "text-right")}
                >
                  <Button
                    type="button"
                    variant="ghost"
                    size="xs"
                    onClick={header.column.getToggleSortingHandler()}
                    className="text-muted-foreground -mx-2 font-medium"
                  >
                    <table.FlexRender header={header} />
                    <SortIcon header={header} />
                  </Button>
                </TableHead>
              ))}
            </TableRow>
          ))}
        </TableHeader>
        <TableBody>
          {table.getRowModel().rows.map((row) => (
            <TableRow key={row.id} className={rowClassName?.(row.original)}>
              {row.getAllCells().map((cell) => (
                <TableCell key={cell.id} className={cn(numeric.has(cell.column.id) && "text-right")}>
                  <table.FlexRender cell={cell} />
                </TableCell>
              ))}
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </div>
  );
}

function SortIcon<Row extends RowData>({ header }: { header: Header<DataFeatures, Row> }) {
  const sorted = header.column.getIsSorted();

  if (sorted === "asc") {
    return <ArrowUp aria-hidden />;
  }

  if (sorted === "desc") {
    return <ArrowDown aria-hidden />;
  }

  return <ArrowUpDown aria-hidden className="opacity-40" />;
}

function ariaSort<Row extends RowData>(header: Header<DataFeatures, Row>): "ascending" | "descending" | "none" {
  const sorted = header.column.getIsSorted();

  if (sorted === "asc") {
    return "ascending";
  }

  return sorted === "desc" ? "descending" : "none";
}
