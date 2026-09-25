import { useCallback } from "react";

import type { DateInterval, YearDate } from "../results/types";
import { truncate } from "../text";
import { PLATE } from "./palette";

export interface DateRow {
  id: string;
  label: string;
  date: YearDate;
  interval: DateInterval | undefined;
  current: boolean;
}

const ROW_HEIGHT = 22;

const LABEL_WIDTH = 170;

const WIDTH = 420;

const AXIS_HEIGHT = 22;

const LABEL_LENGTH = 26;

export function DateIntervals({
  rows,
  onOpen,
}: {
  rows: readonly DateRow[];
  onOpen?: ((id: string) => void) | undefined;
}) {
  const ends = rows
    .flatMap((row) => [row.interval?.lower ?? row.date, row.interval?.upper ?? row.date])
    .toSorted((a, b) => a.year - b.year);

  const first = ends.at(0);
  const last = ends.at(-1);
  const low = first?.year ?? 0;
  const high = last?.year ?? 0;
  const pad = (high - low) * 0.05 || 0.05;
  const scale = { low: low - pad, span: high - low + 2 * pad };
  const height = rows.length * ROW_HEIGHT + AXIS_HEIGHT;

  return (
    <svg viewBox={`0 0 ${WIDTH} ${height}`} width="100%" aria-label="Date in each run">
      <title>Date in each run</title>
      {rows.map((row, index) => (
        <IntervalRow
          key={row.id}
          row={row}
          y={index * ROW_HEIGHT + ROW_HEIGHT / 2}
          scale={scale}
          onOpen={row.current ? undefined : onOpen}
        />
      ))}
      <text x={LABEL_WIDTH} y={height - 6} fontSize={10} fill={PLATE.muted}>
        {first?.date}
      </text>
      <text x={WIDTH - 12} y={height - 6} fontSize={10} fill={PLATE.muted} textAnchor="end">
        {last?.date}
      </text>
    </svg>
  );
}

function IntervalRow({
  row,
  y,
  scale,
  onOpen,
}: {
  row: DateRow;
  y: number;
  scale: { low: number; span: number };
  onOpen: ((id: string) => void) | undefined;
}) {
  const open = useCallback(() => onOpen?.(row.id), [onOpen, row.id]);
  const x = (value: number) => LABEL_WIDTH + ((value - scale.low) / scale.span) * (WIDTH - LABEL_WIDTH - 12);

  return (
    <g
      onClick={onOpen === undefined ? undefined : open}
      className={onOpen === undefined ? undefined : "cursor-pointer"}
    >
      <title>{rowTitle(row)}</title>
      <text
        x={LABEL_WIDTH - 8}
        y={y + 4}
        textAnchor="end"
        fontSize={11}
        fill={row.current ? PLATE.ink : PLATE.muted}
        fontWeight={row.current ? 700 : 400}
      >
        {truncate(row.label, LABEL_LENGTH)}
      </text>
      {row.interval !== undefined && (
        <line
          x1={x(row.interval.lower.year)}
          x2={Math.max(x(row.interval.upper.year), x(row.interval.lower.year) + 0.5)}
          y1={y}
          y2={y}
          stroke={PLATE.accent}
          strokeWidth={3}
          strokeOpacity={0.35}
        />
      )}
      <circle cx={x(row.date.year)} cy={y} r={3.5} fill={row.current ? PLATE.selection : PLATE.accent} />
    </g>
  );
}

function rowTitle(row: DateRow): string {
  if (row.interval === undefined) {
    return `${row.label}: ${row.date.date}`;
  }

  return `${row.label}: ${row.date.date} (${row.interval.lower.date} to ${row.interval.upper.date})`;
}
