import type { RunRecord } from "@neherlab/app-contracts";
import { CircleX, TriangleAlert } from "lucide-react";
import type { ReactNode } from "react";

import { formatDuration } from "../format";
import { Card, CardAction, CardDescription, CardHeader, CardTitle } from "../ui/card";
import { cn } from "../ui/cn";

export function Panel({
  title,
  hint,
  actions,
  children,
  className,
}: {
  title: ReactNode;
  hint?: ReactNode;
  actions?: ReactNode;
  children: ReactNode;
  className?: string;
}) {
  return (
    <Card size="sm" className={cn("min-w-0 gap-0 py-0 shadow-none", className)}>
      <CardHeader className="border-b py-2.5">
        <CardTitle>{title}</CardTitle>
        {hint !== undefined && <CardDescription className="text-xs">{hint}</CardDescription>}
        {actions !== undefined && (
          <CardAction className="flex flex-wrap items-center justify-end gap-1.5">{actions}</CardAction>
        )}
      </CardHeader>
      {children}
    </Card>
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

const TONE_CLASS = {
  caution: "text-warning font-medium",
  fault: "text-destructive font-medium",
} as const;

export function SummaryStrip({ entries }: { entries: readonly SummaryEntry[] }) {
  return (
    <dl className="bg-border grid grid-cols-2 gap-px overflow-hidden rounded-lg border @lg:grid-cols-3 @5xl:grid-cols-6">
      {entries.map((entry) => (
        <div key={entry.label} className="bg-card grid content-start gap-0.5 px-3.5 py-2.5">
          <dt className="text-muted-foreground text-xs">{entry.label}</dt>
          <dd className="text-lg leading-snug font-medium">{entry.value}</dd>
          {entry.detail !== undefined && (
            <dd
              className={cn(
                "flex items-start gap-1 text-xs",
                entry.tone === undefined ? "text-muted-foreground" : TONE_CLASS[entry.tone],
              )}
            >
              {entry.tone === "caution" && <TriangleAlert aria-label="Caution" className="mt-px size-3 shrink-0" />}
              {entry.tone === "fault" && <CircleX aria-label="Fault" className="mt-px size-3 shrink-0" />}
              <span>{entry.detail}</span>
            </dd>
          )}
        </div>
      ))}
    </dl>
  );
}
