import { errorMessage, type RunRecord, type RunResults } from "@neherlab/app-contracts";
import { runsGet } from "@neherlab/app-contracts/client";
import { Alert, AlertDescription, AlertTitle } from "@neherlab/design/ui/alert";
import { Badge } from "@neherlab/design/ui/badge";
import { Button } from "@neherlab/design/ui/button";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@neherlab/design/ui/tabs";
import { Link, Outlet, useChildMatches } from "@tanstack/react-router";
import { useMemo, useState } from "react";
import RotateCcw from "~icons/lucide/rotate-ccw";

import { useRunEvents } from "../api/events";
import { useApi } from "../api/hooks";
import { isNotFound } from "../api/retry";
import { AuspiceStoreSlot, AuspiceStoreSlotContext } from "../auspice/storeSlot";
import { LoadingState, PageShell } from "../components/PageShell";
import { countWarnings, EMPTY_PROGRESS, type RunProgress } from "../results/progress";
import { AncestralResultsView } from "./AncestralResultsView";
import { ClockResultsView } from "./ClockResultsView";
import { LogTab } from "./LogTab";
import { MissingRun } from "./MissingRun";
import { RunContext, useRunOutputs, useRunSummary, useRunView } from "./runContext";
import { RunHeader } from "./RunHeader";
import { RunningView } from "./RunningView";
import { RunWarnings, WarningNames } from "./RunWarnings";
import { SettingsTab } from "./SettingsTab";
import { TimetreeResultsView } from "./TimetreeResultsView";
import { TreeOnlyResults } from "./TreeOnlyResults";
import type { TreeData } from "./TreeView";
import { useRerun } from "./useRerun";
import { VariantDelta } from "./VariantDelta";
import { runViews } from "./views";

const LOG_TAB = "/runs/$id/log";

const LIVE_STATUSES = new Set(["created", "running"]);

export function RunPage({ id }: { id: string }) {
  const { data: record, error } = useApi((context) => runsGet({ ...context, path: { id } }));
  const { data: progress = EMPTY_PROGRESS, error: eventsError } = useRunEvents(id);
  const failure = eventsError === null ? undefined : errorMessage(eventsError);
  const tab = useChildMatches({ select: (matches) => matches.at(0)?.fullPath });
  const [slot] = useState(() => new AuspiceStoreSlot());

  const view = useMemo(
    () => (record === undefined ? undefined : { record, progress, failure }),
    [failure, progress, record],
  );

  if (error !== null && isNotFound(error)) {
    return <MissingRun />;
  }

  if (error !== null) {
    return (
      <PageShell>
        <Alert variant="destructive">
          <AlertTitle>The run cannot be loaded</AlertTitle>
          <AlertDescription>{error.message}</AlertDescription>
        </Alert>
      </PageShell>
    );
  }

  if (record === undefined) {
    return <LoadingState text="Loading the run" />;
  }

  const warnings = countWarnings(progress.entries);

  return (
    <PageShell>
      <RunHeader record={record} />
      <RunWarnings warnings={record.warnings} />
      <Tabs value={tab} className="gap-4">
        <TabsList variant="line" aria-label="Run views" className="border-b">
          {runViews(record.command, record.status).map((entry) => (
            <TabsTrigger
              key={entry.to}
              value={entry.to}
              nativeButton={false}
              render={<Link to={entry.to} params={{ id }} />}
            >
              {entry.label}
              {entry.to === LOG_TAB && warnings > 0 && (
                <Badge variant="outline" className="text-warning border-warning/40">
                  {warnings} {warnings === 1 ? "warning" : "warnings"}
                </Badge>
              )}
            </TabsTrigger>
          ))}
        </TabsList>
        <TabsContent value={tab}>
          <RunContext value={view}>
            <AuspiceStoreSlotContext value={slot}>
              <Outlet />
            </AuspiceStoreSlotContext>
          </RunContext>
        </TabsContent>
      </Tabs>
    </PageShell>
  );
}

export function RunResultsTab() {
  const { record, progress } = useRunView();

  if (LIVE_STATUSES.has(record.status)) {
    return <RunningView record={record} progress={progress} />;
  }

  if (record.status !== "ok") {
    return <EndedRun record={record} progress={progress} />;
  }

  return <FinishedOverview record={record} />;
}

export function RunSettingsTab() {
  const { record } = useRunView();

  const summary = useRunSummary(record.id);

  return <SettingsTab record={record} configName={summary?.config_name} />;
}

export function RunLogTab() {
  const { record, progress, failure } = useRunView();

  const summary = useRunSummary(record.id);

  return <LogTab progress={progress} failure={failure} fileName={summary?.log_name} />;
}

function FinishedOverview({ record }: { record: RunRecord }) {
  const outputs = useRunOutputs(record);

  if (outputs.state === "failed") {
    return (
      <Alert variant="destructive">
        <AlertTitle>The outputs of the run cannot be read</AlertTitle>
        <AlertDescription>{outputs.message}</AlertDescription>
      </Alert>
    );
  }

  if (outputs.state === "loading") {
    return <LoadingState text="Reading the outputs" />;
  }

  const { results, tree } = outputs;

  return (
    <div className="grid gap-4">
      {results.problems.length > 0 && (
        <Alert variant="destructive">
          <AlertTitle>Some outputs cannot be read</AlertTitle>
          <AlertDescription>
            <ul className="list-disc pl-5">
              {results.problems.map((problem) => (
                <li key={problem.path}>
                  <code className="font-mono text-xs">{problem.path}</code>: {problem.message}
                </li>
              ))}
            </ul>
          </AlertDescription>
        </Alert>
      )}
      {results.tree_alignment !== undefined && (
        <Alert>
          <AlertDescription>
            The tree shows the Auspice file of alignment {results.tree_alignment.number} (
            <code className="font-mono text-xs">{results.tree_alignment.path}</code>). The output files hold one Auspice
            file per alignment file.
          </AlertDescription>
        </Alert>
      )}
      <VariantDelta record={record} />
      <CommandResults record={record} results={results} tree={tree} />
    </div>
  );
}

function CommandResults({
  record,
  results,
  tree,
}: {
  record: RunRecord;
  results: RunResults;
  tree: TreeData | undefined;
}) {
  const view = results.results;

  if (view.command === "timetree") {
    return <TimetreeResultsView record={record} data={view.data} tree={tree} />;
  }

  if (view.command === "clock") {
    return <ClockResultsView record={record} data={view.data} tree={tree} />;
  }

  if (view.command === "ancestral") {
    return <AncestralResultsView record={record} data={view.data} tree={tree} />;
  }

  return <TreeOnlyResults record={record} data={view.data} tree={tree} />;
}

function EndedRun({ record, progress }: { record: RunRecord; progress: RunProgress }) {
  const rerun = useRerun(record);
  const error = record.error;
  const affected = record.warnings.filter((warning) => warning.names.length > 0);

  return (
    <div className="grid gap-4">
      <Alert variant="destructive">
        <AlertTitle>{ENDED_TITLES[record.status]}</AlertTitle>
        <AlertDescription className="grid gap-2">
          {error !== undefined && (
            <>
              <p>{error.message}</p>
              {error.causes.length > 0 && (
                <ul className="list-disc pl-5">
                  {error.causes.map((cause) => (
                    <li key={cause}>{cause}</li>
                  ))}
                </ul>
              )}
            </>
          )}
          {affected.map((warning) => (
            <div key={warning.message} className="grid gap-1">
              <p>{warning.message}</p>
              <WarningNames names={warning.names} />
            </div>
          ))}
          <p className="text-muted-foreground">The stages below show the run up to the point where it stopped.</p>
          <Button type="button" variant="outline" size="sm" className="justify-self-start" onClick={rerun}>
            <RotateCcw aria-hidden />
            Edit and run again
          </Button>
        </AlertDescription>
      </Alert>
      <RunningView record={record} progress={progress} />
    </div>
  );
}

const ENDED_TITLES: Readonly<Record<RunRecord["status"], string>> = {
  created: "The run has not started",
  running: "The run is still running",
  ok: "The run finished",
  error: "The run failed",
  cancelled: "The run was cancelled",
  interrupted: "The run was interrupted because the process that ran it stopped",
};
