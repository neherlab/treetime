import { errorMessage } from "@neherlab/app-contracts";
import type { RunRecord } from "@neherlab/app-contracts";
import { runsCancel } from "@neherlab/app-contracts/client";
import { CircleCheck, CircleX } from "lucide-react";
import { type ReactNode, useCallback, useId, useMemo, useState } from "react";

import { useApiMutation } from "../api/hooks";
import { OptionToggle } from "../components/OptionToggle";
import { Panel } from "../components/Panel";
import { formatDuration, formatRate } from "../format";
import { nonFiniteLabel } from "../results/numbers";
import {
  countWarnings,
  filterLog,
  openSections,
  sectionOverrides,
  stageSections,
  type IterationPoint,
  type RunProgress,
  type StageSection,
  type StageState,
} from "../results/progress";
import { Accordion, AccordionContent, AccordionItem, AccordionTrigger } from "../ui/accordion";
import { Badge } from "../ui/badge";
import { Button } from "../ui/button";
import { cn } from "../ui/cn";
import { Progress, ProgressLabel, ProgressValue } from "../ui/progress";
import { Spinner } from "../ui/spinner";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "../ui/table";
import { useToastManager } from "../ui/toast";
import { LogLines } from "./LogLines";
import { RateTrace } from "./RateTrace";

type LiveFilter = "all" | "warnings";

const LIVE_FILTERS: ReadonlyArray<{ value: LiveFilter; label: string }> = [
  { value: "all", label: "All" },
  { value: "warnings", label: "Warnings" },
];

const NO_OVERRIDES: ReadonlyMap<string, boolean> = new Map();

export function RunningView({ record, progress }: { record: RunRecord; progress: RunProgress }) {
  const [filter, setFilter] = useState<LiveFilter>("all");
  const titleId = useId();
  const percent = Math.round(Math.min(1, Math.max(0, progress.fraction)) * 100);

  return (
    <div className="grid gap-4 @3xl:grid-cols-[minmax(0,1fr)_26rem]">
      <section aria-labelledby={titleId} className="grid min-w-0 content-start gap-2">
        <header className="flex flex-wrap items-center gap-x-3 gap-y-1.5">
          <div className="grid gap-0.5">
            <h2 id={titleId} className="font-heading text-sm leading-normal font-bold">
              Log
            </h2>
            <p className="text-muted-foreground text-xs">Open stages follow new lines until you scroll up</p>
          </div>
          <div className="ml-auto">
            <OptionToggle label="Log filter" value={filter} onChange={setFilter} options={LIVE_FILTERS} />
          </div>
        </header>
        <StageLog progress={progress} filter={filter} />
      </section>
      <div className="grid content-start gap-4">
        <Panel title="Progress" hint={progress.message === "" ? undefined : progress.message}>
          <div className="p-3.5">
            <Progress value={percent} aria-label="Run progress">
              <ProgressLabel className="sr-only">Run progress</ProgressLabel>
              <ProgressValue />
            </Progress>
          </div>
        </Panel>
        {progress.iterations.length > 0 && <IterationPanel iterations={progress.iterations} />}
        {record.status === "running" && <CancelButton id={record.id} />}
      </div>
    </div>
  );
}

function StageLog({ progress, filter }: { progress: RunProgress; filter: LiveFilter }) {
  const sections = useMemo(() => stageSections(progress), [progress]);
  const [overrides, setOverrides] = useState<ReadonlyMap<string, boolean>>(NO_OVERRIDES);
  const open = useMemo(() => openSections(sections, overrides), [overrides, sections]);

  const onValueChange = useCallback((next: string[]) => setOverrides(sectionOverrides(sections, next)), [sections]);

  if (sections.length === 0) {
    return <p className="text-muted-foreground py-2 text-xs">No log lines yet.</p>;
  }

  return (
    <Accordion multiple value={open} onValueChange={onValueChange} className="gap-2">
      {sections.map((section) => (
        <StageCard key={section.key} section={section} filter={filter} />
      ))}
    </Accordion>
  );
}

function StageCard({ section, filter }: { section: StageSection; filter: LiveFilter }) {
  const entries = useMemo(() => filterLog(section.entries, filter, ""), [filter, section.entries]);
  const warnings = countWarnings(section.entries);

  return (
    <AccordionItem value={section.key} className="bg-card ring-border rounded-lg ring-1 not-last:border-b-0">
      <AccordionTrigger className="items-center gap-2 px-3 py-2.5 hover:no-underline">
        <StageIcon state={section.state} />
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
          className="max-h-80"
          empty={section.entries.length === 0 ? "No log lines in this stage." : "No matching lines."}
        />
      </AccordionContent>
    </AccordionItem>
  );
}

const STAGE_ICONS = {
  running: <Spinner aria-label="Running" className="text-primary size-3.5 shrink-0" />,
  done: <CircleCheck aria-label="Done" className="text-success size-3.5 shrink-0" />,
  stopped: <CircleX aria-label="Stopped" className="text-destructive size-3.5 shrink-0" />,
} as const satisfies Record<StageState, ReactNode>;

function StageIcon({ state }: { state: StageState }) {
  return STAGE_ICONS[state];
}

function IterationPanel({ iterations }: { iterations: readonly IterationPoint[] }) {
  return (
    <Panel figure title="Clock rate by iteration" hint={latestCaption(iterations.at(-1))}>
      <RateTrace iterations={iterations} />
      <Table className="text-xs">
        <TableHeader>
          <TableRow>
            <TableHead>Iteration</TableHead>
            <TableHead className="text-right">Rate</TableHead>
            <TableHead className="text-right">max Δt</TableHead>
            <TableHead className="text-right">rms Δt</TableHead>
            <TableHead className="text-right">log L</TableHead>
          </TableRow>
        </TableHeader>
        <TableBody>
          {iterations.map((point) => (
            <TableRow key={point.iteration}>
              <TableCell>{point.iteration}</TableCell>
              <NumberCell value={point.clockRate} format={formatRate} />
              <NumberCell value={point.maxTimeChange} format={fixed4} />
              <NumberCell value={point.rmsTimeChange} format={fixed4} />
              <NumberCell value={point.logLhTotal} format={fixed1} />
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </Panel>
  );
}

function NumberCell({ value, format }: { value: number | undefined; format: (value: number) => string }) {
  if (value === undefined) {
    return <TableCell className="text-muted-foreground text-right">-</TableCell>;
  }

  return Number.isFinite(value) ? (
    <TableCell className="text-right">{format(value)}</TableCell>
  ) : (
    <TableCell className="text-destructive text-right font-bold" title="Not a finite number">
      {nonFiniteLabel(value)}
    </TableCell>
  );
}

function CancelButton({ id }: { id: string }) {
  const { mutateAsync: cancelRun } = useApiMutation((context, run: string) =>
    runsCancel({ ...context, path: { id: run } }),
  );

  const toasts = useToastManager();
  const [busy, setBusy] = useState(false);

  const cancel = useCallback(async () => {
    setBusy(true);

    try {
      await cancelRun(id);
    } catch (error: unknown) {
      toasts.add({ title: "The run cannot be cancelled", description: errorMessage(error) });
    } finally {
      setBusy(false);
    }
  }, [cancelRun, id, toasts]);

  const onCancel = useCallback(() => void cancel(), [cancel]);

  return (
    <Button type="button" variant="destructive" onClick={onCancel} disabled={busy}>
      {busy && <Spinner />}
      Cancel run
    </Button>
  );
}

function latestCaption(last: IterationPoint | undefined): string | undefined {
  if (last === undefined) {
    return undefined;
  }

  const rate = `Latest ${formatNumber(last.clockRate, formatRate)}`;

  return last.rSquared === undefined ? rate : `${rate}, R² ${formatNumber(last.rSquared, fixed3)}`;
}

function formatNumber(value: number, format: (value: number) => string): string {
  return Number.isFinite(value) ? format(value) : nonFiniteLabel(value);
}

function fixed4(value: number): string {
  return value.toFixed(4);
}

function fixed3(value: number): string {
  return value.toFixed(3);
}

function fixed1(value: number): string {
  return value.toFixed(1);
}
