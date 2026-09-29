import { type ReactNode, useCallback, useId, useMemo, useState } from "react";
import CircleCheck from "~icons/lucide/circle-check";
import CircleX from "~icons/lucide/circle-x";

import { formatDuration } from "../format";
import {
  countWarnings,
  filterLog,
  openSections,
  sectionOverrides,
  stageSections,
  type LogFilter,
  type RunProgress,
  type StageSection,
  type StageState,
} from "../results/progress";
import { Accordion, AccordionContent, AccordionItem, AccordionTrigger } from "../ui/accordion";
import { Badge } from "../ui/badge";
import { cn } from "../ui/cn";
import { Spinner } from "../ui/spinner";
import { LogLines } from "./LogLines";

const NO_OVERRIDES: ReadonlyMap<string, boolean> = new Map();

const STAGE_ICONS = {
  running: <Spinner aria-label="Running" className="text-primary size-3.5 shrink-0" />,
  done: <CircleCheck aria-label="Done" className="text-success size-3.5 shrink-0" />,
  stopped: <CircleX aria-label="Stopped" className="text-destructive size-3.5 shrink-0" />,
} as const satisfies Record<StageState, ReactNode>;

export function StageLog({
  progress,
  filter,
  query,
  hint,
  actions,
  children,
}: {
  progress: RunProgress;
  filter: LogFilter;
  query: string;
  hint: string;
  actions: ReactNode;
  children?: ReactNode;
}) {
  const titleId = useId();
  const sections = useMemo(() => stageSections(progress), [progress]);
  const [overrides, setOverrides] = useState<ReadonlyMap<string, boolean>>(NO_OVERRIDES);
  const open = useMemo(() => openSections(sections, overrides), [overrides, sections]);

  const onValueChange = useCallback((next: string[]) => setOverrides(sectionOverrides(sections, next)), [sections]);

  return (
    <section aria-labelledby={titleId} className="grid min-w-0 content-start gap-2">
      <header className="flex flex-wrap items-center gap-x-3 gap-y-1.5">
        <div className="grid gap-0.5">
          <h2 id={titleId} className="font-heading text-sm leading-normal font-bold">
            Log
          </h2>
          <p className="text-muted-foreground text-xs">{hint}</p>
        </div>
        <div className="ml-auto flex flex-wrap items-center justify-end gap-1.5">{actions}</div>
      </header>
      {children}
      {sections.length === 0 ? (
        <p className="text-muted-foreground py-2 text-xs">No log lines yet.</p>
      ) : (
        <Accordion multiple value={open} onValueChange={onValueChange} className="gap-2">
          {sections.map((section) => (
            <StageCard
              key={section.key}
              section={section}
              filter={filter}
              query={query}
              follow={progress.terminal === undefined}
            />
          ))}
        </Accordion>
      )}
    </section>
  );
}

function StageCard({
  section,
  filter,
  query,
  follow,
}: {
  section: StageSection;
  filter: LogFilter;
  query: string;
  follow: boolean;
}) {
  const entries = useMemo(() => filterLog(section.entries, filter, query), [filter, query, section.entries]);
  const warnings = countWarnings(section.entries);

  return (
    <AccordionItem value={section.key} className="bg-card ring-border rounded-lg ring-1 not-last:border-b-0">
      <AccordionTrigger className="items-center gap-2 px-3 py-2.5 hover:no-underline **:data-[slot=accordion-trigger-icon]:ml-0">
        {STAGE_ICONS[section.state]}
        <span className={cn("min-w-0 truncate", section.state !== "running" && "font-normal")}>{section.name}</span>
        <span className="text-muted-foreground ml-auto flex shrink-0 items-center gap-2 text-xs font-normal">
          {warnings > 0 && (
            <Badge variant="outline" className="text-warning">
              {warnings === 1 ? "1 warning" : `${warnings} warnings`}
            </Badge>
          )}
          <span>{section.entries.length === 1 ? "1 line" : `${section.entries.length} lines`}</span>
          <span className="w-14 text-right">{formatDuration(section.seconds)}</span>
        </span>
      </AccordionTrigger>
      <AccordionContent className="px-2.5 pb-2.5">
        <LogLines
          entries={entries}
          follow={follow}
          className="max-h-80"
          empty={section.entries.length === 0 ? "No log lines in this stage." : "No matching lines."}
        />
      </AccordionContent>
    </AccordionItem>
  );
}
