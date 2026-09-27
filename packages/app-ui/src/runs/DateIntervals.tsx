import { useCallback, useMemo, useState } from "react";

import { useElementWidth } from "../hooks/useElementWidth";
import type { DateInterval, YearDate } from "../results/types";
import { cn } from "../ui/cn";
import { CHART, clamp, niceAxis, yearTick } from "./palette";

export interface DateRow {
  id: string;
  label: string;
  date: YearDate;
  interval: DateInterval | undefined;
  current: boolean;
}

const TRACK_HEIGHT = 22;

const AXIS_HEIGHT = 20;

const PX_PER_TICK = 80;

const MIN_TICKS = 2;

const MAX_TICKS = 6;

const MIN_SPAN_YEARS = 0.05;

export function DateIntervals({
  rows,
  onOpen,
}: {
  rows: readonly DateRow[];
  onOpen?: ((id: string) => void) | undefined;
}) {
  const [axis, setAxis] = useState<HTMLDivElement | null>(null);
  const tickCount = clamp(Math.floor(useElementWidth(axis) / PX_PER_TICK), MIN_TICKS, MAX_TICKS);
  const scale = useMemo(() => dateScale(rows, tickCount), [rows, tickCount]);

  return (
    <div className="grid grid-cols-[minmax(0,12rem)_auto_minmax(0,1fr)] items-center gap-x-3 text-xs">
      {rows.map((row) => (
        <IntervalRow key={row.id} row={row} scale={scale} onOpen={row.current ? undefined : onOpen} />
      ))}
      <span className="col-span-2" />
      <div ref={setAxis}>
        <svg width="100%" height={AXIS_HEIGHT} className="overflow-visible" aria-hidden>
          <line x1="0%" x2="100%" y1={0.5} y2={0.5} stroke={CHART.faint} />
          {scale.ticks.map((tick, index) => (
            <g key={tick}>
              <line x1={scale.at(tick)} x2={scale.at(tick)} y1={0} y2={4} stroke={CHART.faint} />
              <text
                x={scale.at(tick)}
                y={AXIS_HEIGHT - 3}
                fontSize={10}
                fill={CHART.muted}
                textAnchor={tickAnchor(index, scale.ticks.length)}
              >
                {yearTick(tick)}
              </text>
            </g>
          ))}
        </svg>
      </div>
    </div>
  );
}

function IntervalRow({
  row,
  scale,
  onOpen,
}: {
  row: DateRow;
  scale: DateScale;
  onOpen: ((id: string) => void) | undefined;
}) {
  const open = useCallback(() => onOpen?.(row.id), [onOpen, row.id]);
  const title = rowTitle(row);
  const labelClass = cn("truncate text-left", row.current ? "text-ink font-bold" : "text-ink-muted");

  return (
    <>
      {onOpen === undefined ? (
        <span className={labelClass} title={title}>
          {row.label}
        </span>
      ) : (
        <button
          type="button"
          className={cn(labelClass, "hover:text-accent cursor-pointer underline-offset-2 hover:underline")}
          title={`${title}. Open this run`}
          onClick={open}
        >
          {row.label}
        </button>
      )}
      <span className={cn("tabular-nums", row.current ? "text-ink font-bold" : "text-ink-muted")}>{row.date.date}</span>
      <svg width="100%" height={TRACK_HEIGHT} className="overflow-visible" aria-hidden>
        <title>{title}</title>
        <line x1="0%" x2="100%" y1={TRACK_HEIGHT / 2} y2={TRACK_HEIGHT / 2} stroke={CHART.grid} />
        {row.interval !== undefined && (
          <line
            x1={scale.at(row.interval.lower.year)}
            x2={scale.at(row.interval.upper.year)}
            y1={TRACK_HEIGHT / 2}
            y2={TRACK_HEIGHT / 2}
            stroke={CHART.accent}
            strokeWidth={4}
            strokeOpacity={0.35}
            strokeLinecap="round"
          />
        )}
        <circle
          cx={scale.at(row.date.year)}
          cy={TRACK_HEIGHT / 2}
          r={4}
          fill={row.current ? CHART.selection : CHART.accent}
        />
      </svg>
    </>
  );
}

function dateScale(rows: readonly DateRow[], tickCount: number): DateScale {
  const ends = rows.flatMap((row) => [
    row.interval?.lower.year ?? row.date.year,
    row.interval?.upper.year ?? row.date.year,
  ]);

  const low = Math.min(...ends);
  const high = Math.max(...ends);
  const pad = Math.max((high - low) * 0.05, MIN_SPAN_YEARS / 2);

  const {
    domain: [from, to],
    ticks,
  } = niceAxis([low - pad, high + pad], tickCount);

  return { ticks, at: (year: number) => `${((year - from) / (to - from)) * 100}%` };
}

function tickAnchor(index: number, count: number): "start" | "middle" | "end" {
  if (index === 0) {
    return "start";
  }

  return index === count - 1 ? "end" : "middle";
}

interface DateScale {
  ticks: number[];
  at: (year: number) => string;
}

function rowTitle(row: DateRow): string {
  if (row.interval === undefined) {
    return `${row.label}: ${row.date.date}`;
  }

  return `${row.label}: ${row.date.date} (${row.interval.lower.date} to ${row.interval.upper.date})`;
}
