import { errorMessage } from "@neherlab/app-contracts";
import type { RunRecord, RunResults } from "@neherlab/app-contracts";
import { runsAuspice, runsGet, runsList, runsResults, runsUpdate } from "@neherlab/app-contracts/client";
import { Link, useNavigate } from "@tanstack/react-router";
import { Pin, PinOff, RotateCcw } from "lucide-react";
import { DateTime } from "luxon";
import { useCallback, useMemo } from "react";

import { useRunEvents } from "../api/events";
import { useApi, useApiMutation } from "../api/hooks";
import { LoadingState, PageShell } from "../components/PageShell";
import { formatDuration } from "../format";
import { countWarnings, EMPTY_PROGRESS, type RunProgress } from "../results/progress";
import { COMMAND_INFO } from "../settings/commands";
import { StatusIcon, statusLabel } from "../shell/StatusIcon";
import { Alert, AlertDescription, AlertTitle } from "../ui/alert";
import { Badge } from "../ui/badge";
import { Button } from "../ui/button";
import { NativeSelect, NativeSelectOption } from "../ui/native-select";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "../ui/tabs";
import { useToastManager } from "../ui/toast";
import { AncestralResults } from "./AncestralResults";
import { ClockResults } from "./ClockResults";
import { LogTab } from "./LogTab";
import { MugrationResults } from "./MugrationResults";
import { RunningView } from "./RunningView";
import { RunTitle } from "./RunTitle";
import { SettingsTab } from "./SettingsTab";
import { TimetreeResults } from "./TimetreeResults";
import { TreeOnlyResults } from "./TreeOnlyResults";
import type { TreeData } from "./TreeView";
import { useRerun } from "./useRerun";

type RunTab = "results" | "settings" | "log";

const TABS: ReadonlyArray<{
  tab: RunTab;
  label: string;
  to: "/runs/$id/results" | "/runs/$id/settings" | "/runs/$id/log";
}> = [
  { tab: "results", label: "Results", to: "/runs/$id/results" },
  { tab: "settings", label: "Settings", to: "/runs/$id/settings" },
  { tab: "log", label: "Log", to: "/runs/$id/log" },
];

const LIVE_STATUSES = new Set(["created", "running"]);

export function RunPage({ id, tab }: { id: string; tab: RunTab }) {
  const { data: record, error } = useApi((context) => runsGet({ ...context, path: { id } }));
  const { data: progress = EMPTY_PROGRESS, error: eventsError } = useRunEvents(id);
  const failure = eventsError === null ? undefined : errorMessage(eventsError);

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
      <Tabs value={tab} className="gap-4">
        <TabsList variant="line" aria-label="Run views" className="border-b">
          {TABS.map((entry) => (
            <TabsTrigger
              key={entry.tab}
              value={entry.tab}
              nativeButton={false}
              render={<Link to={entry.to} params={{ id }} />}
            >
              {entry.label}
              {entry.tab === "log" && warnings > 0 && (
                <Badge variant="outline" className="text-warning border-warning/40">
                  {warnings} {warnings === 1 ? "warning" : "warnings"}
                </Badge>
              )}
            </TabsTrigger>
          ))}
        </TabsList>
        <TabsContent value={tab}>
          {tab === "results" && <ResultsTab record={record} progress={progress} />}
          {tab === "settings" && <SettingsTab record={record} />}
          {tab === "log" && <LogTab progress={progress} failure={failure} />}
        </TabsContent>
      </Tabs>
    </PageShell>
  );
}

function ResultsTab({ record, progress }: { record: RunRecord; progress: RunProgress }) {
  if (LIVE_STATUSES.has(record.status)) {
    return <RunningView record={record} progress={progress} />;
  }

  if (record.status !== "ok") {
    return <EndedRun record={record} progress={progress} />;
  }

  return <FinishedResults record={record} />;
}

function FinishedResults({ record }: { record: RunRecord }) {
  const { data: results, error } = useApi((context) => runsResults({ ...context, path: { id: record.id } }), {
    staleTime: Infinity,
  });

  const hasTree = results?.tree !== null && results?.tree !== undefined;

  const { data: document, error: treeError } = useApi(
    (context) => runsAuspice({ ...context, path: { id: record.id } }),
    { enabled: hasTree, staleTime: Infinity },
  );

  const tree = useMemo<TreeData | undefined>(
    () =>
      results?.tree === null || results?.tree === undefined || document === undefined
        ? undefined
        : { document, tree: results.tree },
    [document, results],
  );

  const failure = error ?? treeError;

  if (failure !== null) {
    return (
      <Alert variant="destructive">
        <AlertTitle>The outputs of the run cannot be read</AlertTitle>
        <AlertDescription>{failure.message}</AlertDescription>
      </Alert>
    );
  }

  if (results === undefined || (hasTree && document === undefined)) {
    return <LoadingState text="Reading the outputs" />;
  }

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
    return <TimetreeResults record={record} results={results} data={view.data} tree={tree} />;
  }

  if (view.command === "clock") {
    return <ClockResults record={record} results={results} data={view.data} tree={tree} />;
  }

  if (view.command === "ancestral") {
    return <AncestralResults record={record} results={results} data={view.data} tree={tree} />;
  }

  if (view.command === "mugration") {
    return <MugrationResults record={record} results={results} data={view.data} tree={tree} />;
  }

  return <TreeOnlyResults record={record} results={results} data={view.data} tree={tree} />;
}

function EndedRun({ record, progress }: { record: RunRecord; progress: RunProgress }) {
  const rerun = useRerun(record);
  const error = record.error;

  return (
    <div className="grid gap-4">
      <Alert variant="destructive">
        <AlertTitle>{ENDED_TITLES[record.status]}</AlertTitle>
        <AlertDescription className="grid gap-2">
          {error !== null && error !== undefined && (
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
          <p className="text-muted-foreground">
            The progress and log below show the run up to the point where it stopped.
          </p>
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

function RunHeader({ record }: { record: RunRecord }) {
  const navigate = useNavigate();
  const toasts = useToastManager();
  const rerun = useRerun(record);
  const { data: list } = useApi((context) => runsList(context));

  const { mutateAsync: updateRun } = useApiMutation((context, pinned: boolean) =>
    runsUpdate({ ...context, path: { id: record.id }, body: { pinned } }),
  );

  const created = DateTime.fromISO(record.created_at);

  const comparable = useMemo(
    () => (list?.runs ?? []).filter((run) => run.id !== record.id && run.status === "ok"),
    [list, record.id],
  );

  const togglePin = useCallback(async () => {
    try {
      await updateRun(!record.pinned);
    } catch (error: unknown) {
      toasts.add({ title: "The run cannot be pinned", description: errorMessage(error) });
    }
  }, [record.pinned, toasts, updateRun]);

  const onTogglePin = useCallback(() => void togglePin(), [togglePin]);

  const onCompare = useCallback(
    (event: React.ChangeEvent<HTMLSelectElement>) => {
      const other = event.target.value;

      if (other !== "") {
        void navigate({ to: "/compare/$a/$b", params: { a: record.id, b: other } });
      }
    },
    [navigate, record.id],
  );

  return (
    <header className="flex flex-wrap items-start gap-4">
      <div className="grid min-w-0 gap-1">
        <RunTitle record={record} />
        <p className="text-muted-foreground flex flex-wrap items-center gap-x-3 gap-y-1 text-sm">
          <span className="text-foreground inline-flex items-center gap-1.5 font-medium">
            <StatusIcon status={record.status} />
            {statusLabel(record.status)}
          </span>
          <span>{COMMAND_INFO[record.command].label}</span>
          <span>{created.isValid ? created.toFormat("d LLL yyyy, HH:mm") : record.created_at}</span>
          {record.duration_seconds !== null && record.duration_seconds !== undefined && (
            <span>{formatDuration(record.duration_seconds)}</span>
          )}
          <span>TreeTime {record.treetime_version}</span>
        </p>
      </div>
      <div className="ml-auto flex flex-wrap items-center gap-1.5">
        {comparable.length > 0 && (
          <NativeSelect size="sm" aria-label="Compare with another run" value="" onChange={onCompare} className="w-44">
            <NativeSelectOption value="">Compare with...</NativeSelectOption>
            {comparable.map((run) => (
              <NativeSelectOption key={run.id} value={run.id}>
                {run.title}
              </NativeSelectOption>
            ))}
          </NativeSelect>
        )}
        <Button type="button" variant="ghost" size="sm" onClick={onTogglePin}>
          {record.pinned ? <PinOff aria-hidden /> : <Pin aria-hidden />}
          {record.pinned ? "Unpin" : "Pin"}
        </Button>
        <Button type="button" variant="outline" size="sm" onClick={rerun}>
          <RotateCcw aria-hidden />
          Edit and run again
        </Button>
      </div>
    </header>
  );
}
