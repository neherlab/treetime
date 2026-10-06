import type { RecurrentMutation } from "@neherlab/app-contracts";
import { useCallback, useMemo } from "react";

import { DataTable, dataColumns } from "../components/DataTable";
import { Button } from "../ui/button";
import { drmText } from "./homoplasy";

const column = dataColumns<RecurrentMutation>();

const NUMERIC = new Set(["branches", "terminal"]);

const BRANCHES_SORT = [{ id: "branches", desc: true }];

export function RecurrentTable({
  label,
  rows,
  pressed,
  drmAnnotated,
  onColor,
}: {
  label: string;
  rows: readonly RecurrentMutation[];
  pressed: ReadonlySet<string>;
  drmAnnotated: boolean;
  onColor: (position: number) => void;
}) {
  const columns = useMemo(() => recurrentColumns(pressed, drmAnnotated, onColor), [drmAnnotated, onColor, pressed]);

  const rowClassName = useCallback(
    (row: RecurrentMutation) => (pressed.has(row.mutation) ? "bg-accent" : undefined),
    [pressed],
  );

  return (
    <DataTable
      label={label}
      columns={columns}
      rows={rows}
      rowId={mutationKey}
      initialSorting={BRANCHES_SORT}
      numeric={NUMERIC}
      rowClassName={rowClassName}
    />
  );
}

export function PositionButton({
  position,
  text,
  label,
  pressed,
  onColor,
}: {
  position: number;
  text: string;
  label: string;
  pressed: boolean;
  onColor: (position: number) => void;
}) {
  const color = useCallback(() => onColor(position), [onColor, position]);

  return (
    <Button
      type="button"
      variant="ghost"
      size="xs"
      className="font-mono"
      aria-label={label}
      aria-pressed={pressed}
      onClick={color}
    >
      {text}
    </Button>
  );
}

function recurrentColumns(pressed: ReadonlySet<string>, drmAnnotated: boolean, onColor: (position: number) => void) {
  return column.columns([
    column.accessor((row) => row.mutation, {
      id: "mutation",
      header: "Mutation",
      cell: ({ row }) => (
        <PositionButton
          position={row.original.position}
          text={row.original.mutation}
          label={`Color the tree by ${row.original.mutation}`}
          pressed={pressed.has(row.original.mutation)}
          onColor={onColor}
        />
      ),
    }),
    column.accessor((row) => row.branches, { id: "branches", header: "Branches" }),
    column.accessor((row) => row.terminal_branches, { id: "terminal", header: "Terminal" }),
    ...(drmAnnotated
      ? [column.accessor((row) => (row.drm === undefined ? "" : drmText(row.drm)), { id: "drm", header: "DRM" })]
      : []),
  ]);
}

function mutationKey(row: RecurrentMutation): string {
  return row.mutation;
}
