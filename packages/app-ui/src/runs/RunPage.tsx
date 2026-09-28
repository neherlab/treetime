import { errorMessage } from "@neherlab/app-contracts";
import type { RunRecord, RunResults } from "@neherlab/app-contracts";
import { runsAuspice, runsGet, runsList, runsResults, runsUpdate } from "@neherlab/app-contracts/client";
import { Link, useNavigate } from "@tanstack/react-router";
import { DateTime } from "luxon";
import { useCallback, useMemo } from "react";

import { useRunEvents } from "../api/events";
import { useApi, useApiMutation } from "../api/hooks";
import { formatDuration } from "../format";
import { countWarnings, EMPTY_PROGRESS, type RunProgress } from "../results/progress";
import { COMMAND_INFO } from "../settings/commands";
import { StatusIcon, statusLabel } from "../shell/StatusIcon";
import { Button, Toast } from "../ui";
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
    return <p className="text-signal-danger p-10 text-center">The run cannot be loaded: {error.message}</p>;
  }

  if (record === undefined) {
    return <p className="text-ink-muted p-10 text-center">Loading the run...</p>;
  }

  const warnings = countWarnings(progress.entries);

  return (
    <div className="mx-auto max-w-[110rem] px-5 pt-4 pb-16">
      <RunHeader record={record} />
      <nav aria-label="Run views" className="border-line mb-4 flex gap-0.5 border-b">
        {TABS.map((entry) => (
          <Link
            key={entry.tab}
            to={entry.to}
            params={{ id }}
            aria-current={entry.tab === tab ? "page" : undefined}
            className="text-ink-muted aria-[current=page]:border-accent aria-[current=page]:text-ink -mb-px border-b-2 border-transparent px-3 py-2 font-bold"
          >
            {entry.label}
            {entry.tab === "log" && warnings > 0 && (
              <span className="bg-signal-warn-subtle text-signal-warn ml-1.5 rounded-sm px-1.5 py-0.5 text-xs">
                {warnings} {warnings === 1 ? "warning" : "warnings"}
              </span>
            )}
          </Link>
        ))}
      </nav>
      {tab === "results" && <ResultsTab record={record} progress={progress} />}
      {tab === "settings" && <SettingsTab record={record} />}
      {tab === "log" && <LogTab progress={progress} failure={failure} />}
    </div>
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

  if (error !== null || treeError !== null) {
    return <p className="text-signal-danger">The outputs of the run cannot be read: {(error ?? treeError)?.message}</p>;
  }

  if (results === undefined || (hasTree && document === undefined)) {
    return <p className="text-ink-muted p-10 text-center">Reading the outputs...</p>;
  }

  return (
    <div className="grid gap-3.5">
      {results.problems.length > 0 && (
        <div role="alert" className="border-signal-danger bg-signal-danger-subtle rounded-lg border px-4 py-3">
          <h3 className="mb-1 font-bold">Some outputs cannot be read</h3>
          <ul className="m-0 list-disc pl-5">
            {results.problems.map((problem) => (
              <li key={problem.path}>
                <code className="font-mono text-xs">{problem.path}</code>: {problem.message}
              </li>
            ))}
          </ul>
        </div>
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
    <div className="grid gap-3.5">
      <div role="alert" className="border-signal-danger bg-signal-danger-subtle rounded-lg border px-4 py-3.5">
        <h3 className="mb-1.5 text-[0.9375rem] font-bold">{endedTitle(record)}</h3>
        {error !== null && error !== undefined && (
          <>
            <p className="m-0">{error.message}</p>
            {error.causes.length > 0 && (
              <ul className="text-ink-muted mt-1 list-disc pl-5">
                {error.causes.map((cause) => (
                  <li key={cause}>{cause}</li>
                ))}
              </ul>
            )}
          </>
        )}
        <p className="text-ink-muted mt-2 mb-0">
          The progress and log below show the run up to the point where it stopped.
        </p>
        <Button type="button" variant="outline" size="sm" className="mt-2.5" onClick={rerun}>
          Edit and run again
        </Button>
      </div>
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

function endedTitle(record: RunRecord): string {
  return ENDED_TITLES[record.status];
}

function RunHeader({ record }: { record: RunRecord }) {
  const navigate = useNavigate();
  const toasts = Toast.useToastManager();
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
    <div className="mb-3.5 flex flex-wrap items-start gap-4">
      <div className="min-w-0">
        <RunTitle record={record} />
        <div className="text-ink-muted mt-1 flex flex-wrap items-center gap-x-3 gap-y-1.5">
          <span className="inline-flex items-center gap-1.5 font-bold">
            <StatusIcon status={record.status} />
            {statusLabel(record.status)}
          </span>
          <span>{COMMAND_INFO[record.command].label}</span>
          <span>{created.isValid ? created.toFormat("d LLL yyyy, HH:mm") : record.created_at}</span>
          {record.duration_seconds !== null && record.duration_seconds !== undefined && (
            <span>{formatDuration(record.duration_seconds)}</span>
          )}
          <span>TreeTime {record.treetime_version}</span>
        </div>
      </div>
      <div className="ml-auto flex flex-wrap gap-1.5">
        {comparable.length > 0 && (
          <select
            aria-label="Compare with another run"
            value=""
            onChange={onCompare}
            className="border-line-strong bg-surface-1 h-7 rounded-md border px-2 text-xs"
          >
            <option value="">Compare with...</option>
            {comparable.map((run) => (
              <option key={run.id} value={run.id}>
                {run.title}
              </option>
            ))}
          </select>
        )}
        <Button type="button" variant="ghost" size="sm" onClick={onTogglePin}>
          {record.pinned ? "Unpin" : "Pin"}
        </Button>
        <Button type="button" variant="outline" size="sm" onClick={rerun}>
          Edit and run again
        </Button>
      </div>
    </div>
  );
}
