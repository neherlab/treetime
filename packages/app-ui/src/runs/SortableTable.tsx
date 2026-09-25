import { ArrowDown, ArrowUp } from "lucide-react";
import { useCallback, useMemo, useState, type ReactNode } from "react";

import { cn } from "../ui";

export type Column<Row> = NumberColumn<Row> | TextColumn<Row>;

interface NumberColumn<Row> {
  kind: "number";
  key: string;
  label: string;
  value: (row: Row) => number;
  render?: (row: Row) => ReactNode;
}

interface TextColumn<Row> {
  kind: "text";
  key: string;
  label: string;
  value: (row: Row) => string;
  render?: (row: Row) => ReactNode;
}

interface Sort {
  key: string;
  descending: boolean;
}

export function SortableTable<Row>({
  columns,
  rows,
  rowKey,
  initialSort,
  rowTone,
  label,
}: {
  columns: ReadonlyArray<Column<Row>>;
  rows: readonly Row[];
  rowKey: (row: Row) => string;
  initialSort: Sort;
  rowTone?: ((row: Row) => "caution" | undefined) | undefined;
  label: string;
}) {
  const [sort, setSort] = useState<Sort>(initialSort);
  const column = columns.find((candidate) => candidate.key === sort.key);

  const sorted = useMemo(() => {
    if (column === undefined) {
      return rows;
    }

    const direction = sort.descending ? -1 : 1;

    return rows.toSorted((left, right) => compareRows(column, left, right) * direction);
  }, [column, rows, sort.descending]);

  return (
    <div className="max-h-[36rem] overflow-auto">
      <table className="w-full border-collapse text-left text-xs" aria-label={label}>
        <thead className="bg-surface-1 sticky top-0">
          <tr>
            {columns.map((entry) => (
              <HeaderCell key={entry.key} column={entry} sort={sort} onSort={setSort} />
            ))}
          </tr>
        </thead>
        <tbody>
          {sorted.map((row) => (
            <tr
              key={rowKey(row)}
              className={cn("border-line border-t", rowTone?.(row) === "caution" && "bg-signal-warn-subtle")}
            >
              {columns.map((entry) => (
                <td key={entry.key} className={cn("px-3 py-1", entry.kind === "number" && "text-right tabular-nums")}>
                  {entry.render?.(row) ?? String(entry.value(row))}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

function HeaderCell<Row>({ column, sort, onSort }: { column: Column<Row>; sort: Sort; onSort: (sort: Sort) => void }) {
  const active = sort.key === column.key;

  const toggle = useCallback(
    () => onSort({ key: column.key, descending: active ? !sort.descending : column.kind === "number" }),
    [active, column.key, column.kind, onSort, sort.descending],
  );

  return (
    <th
      className={cn("text-ink-faint px-3 py-1.5 font-normal", column.kind === "number" && "text-right")}
      aria-sort={active ? (sort.descending ? "descending" : "ascending") : "none"}
    >
      <button type="button" onClick={toggle} className="hover:text-ink inline-flex items-center gap-1">
        {column.label}
        {active && (sort.descending ? <ArrowDown size={11} aria-hidden /> : <ArrowUp size={11} aria-hidden />)}
      </button>
    </th>
  );
}

function compareRows<Row>(column: Column<Row>, left: Row, right: Row): number {
  return column.kind === "number"
    ? column.value(left) - column.value(right)
    : column.value(left).localeCompare(column.value(right));
}
