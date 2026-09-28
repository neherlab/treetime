import type { RunRecord } from "@neherlab/app-contracts";
import { CircleX, TriangleAlert } from "lucide-react";
import type { ReactNode } from "react";

import { formatDuration } from "../format";
import { cn } from "../ui";

export function Panel({
  title,
  hint,
  actions,
  children,
  className,
}: {
  title: string;
  hint?: ReactNode;
  actions?: ReactNode;
  children: ReactNode;
  className?: string;
}) {
  return (
    <section className={cn("border-line bg-surface-1 min-w-0 rounded-lg border", className)}>
      <header className="border-line flex flex-wrap items-center gap-x-2.5 gap-y-1 border-b px-3.5 py-2.5">
        <h3 className="font-bold">{title}</h3>
        {hint !== undefined && <span className="text-ink-faint text-xs">{hint}</span>}
        {actions !== undefined && <div className="ml-auto flex flex-wrap items-center gap-1.5">{actions}</div>}
      </header>
      {children}
    </section>
  );
}

export interface SummaryEntry {
  label: string;
  value: ReactNode;
  detail?: ReactNode;
  tone?: "caution" | "fault" | undefined;
}

export function runTimeEntry(record: RunRecord): SummaryEntry {
  return {
    label: "Run time",
    value:
      record.duration_seconds === null || record.duration_seconds === undefined
        ? "-"
        : formatDuration(record.duration_seconds),
  };
}

export function SummaryStrip({ entries }: { entries: readonly SummaryEntry[] }) {
  return (
    <dl className="border-line bg-surface-1 m-0 grid grid-cols-2 gap-px overflow-hidden rounded-lg border @lg:grid-cols-3 @5xl:grid-cols-6">
      {entries.map((entry) => (
        <div key={entry.label} className="bg-surface-1 px-3.5 py-2.5">
          <dt className="text-ink-faint text-xs">{entry.label}</dt>
          <dd className="m-0 text-lg leading-snug">{entry.value}</dd>
          {entry.detail !== undefined && (
            <dd
              className={cn(
                "m-0 text-xs",
                entry.tone === undefined && "text-ink-muted",
                entry.tone === "caution" && "text-signal-warn font-bold",
                entry.tone === "fault" && "text-signal-danger font-bold",
              )}
            >
              {entry.tone === "caution" && (
                <TriangleAlert size={12} aria-label="Caution" className="mr-1 inline align-[-2px]" />
              )}
              {entry.tone === "fault" && <CircleX size={12} aria-label="Fault" className="mr-1 inline align-[-2px]" />}
              {entry.detail}
            </dd>
          )}
        </div>
      ))}
    </dl>
  );
}
