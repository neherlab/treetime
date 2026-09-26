import type { AppCommand, RunSummaryResult } from "@neherlab/app-contracts";
import { runsList } from "@neherlab/app-contracts/client";
import { Link, useNavigate } from "@tanstack/react-router";
import { Plus } from "lucide-react";
import { DateTime } from "luxon";
import { useCallback, useMemo } from "react";

import { useApi } from "../api/hooks";
import { headlineText } from "../format";
import { APP_COMMANDS, COMMAND_INFO } from "../settings/commands";
import { useShellStore } from "../store/shell";
import { truncate } from "../text";
import { Button, cn } from "../ui";
import { changedFlags, groupRuns, listedRuns } from "./runList";
import { StatusIcon, statusLabel } from "./StatusIcon";
import { useCurrentRunId } from "./useCurrentRunId";

export const RUN_FILTER_ID = "run-filter";

const CHIPS_SHOWN = 2;

const CHIP_LENGTH = 26;

export function Sidebar() {
  const { data, error } = useApi((context) => runsList(context));
  const runFilter = useShellStore((state) => state.runFilter);
  const setRunFilter = useShellStore((state) => state.setRunFilter);
  const commandFilter = useShellStore((state) => state.commandFilter);
  const navigate = useNavigate();
  const allRuns = useMemo(() => listedRuns(data?.runs ?? [], "", null), [data]);
  const runs = useMemo(() => listedRuns(allRuns, runFilter, commandFilter), [allRuns, runFilter, commandFilter]);
  const groups = useMemo(() => groupRuns(runs, DateTime.now()), [runs]);

  const commands = useMemo(
    () => APP_COMMANDS.filter((command) => allRuns.some((run) => run.command === command)),
    [allRuns],
  );

  const openNew = useCallback(() => void navigate({ to: "/new" }), [navigate]);

  const onFilter = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => setRunFilter(event.target.value),
    [setRunFilter],
  );

  return (
    <nav aria-label="Runs" className="border-line bg-surface-1 hidden min-h-0 flex-col border-r md:flex">
      <div className="grid gap-2 px-3 pt-3 pb-2">
        <Button className="w-full" onClick={openNew}>
          <Plus size={14} aria-hidden />
          New analysis
          <kbd className="bg-accent-fg/20 rounded-sm px-1 font-mono text-[0.6875rem]">N</kbd>
        </Button>
        <input
          id={RUN_FILTER_ID}
          type="search"
          value={runFilter}
          onChange={onFilter}
          placeholder="Filter runs by name or setting"
          aria-label="Filter runs"
          className="border-line bg-surface-2 focus-visible:ring-accent rounded-md border px-2.5 py-1.5 outline-none focus-visible:ring-2"
        />
        {commands.length > 1 && <CommandFilter commands={commands} />}
      </div>
      <div className="min-h-0 flex-1 overflow-auto px-1.5 pb-3">
        {error !== null && <p className="text-signal-danger p-3 text-xs">The runs cannot be listed: {error.message}</p>}
        {groups.map((group) => (
          <section key={group.label} aria-label={group.label}>
            <h2 className="text-ink-faint flex justify-between px-2 pt-3 pb-1 text-xs font-bold">
              <span>{group.label}</span>
              <span>{group.runs.length}</span>
            </h2>
            {group.runs.map((run) => (
              <RunRow key={run.id} run={run} />
            ))}
          </section>
        ))}
        {groups.length === 0 && data !== undefined && (
          <p className="text-ink-muted p-6 text-center">
            {allRuns.length === 0 ? "No runs yet. Start a new analysis." : "No run matches the filter."}
          </p>
        )}
      </div>
      <CompareBar />
    </nav>
  );
}

function CommandFilter({ commands }: { commands: readonly AppCommand[] }) {
  return (
    <fieldset className="m-0 flex flex-wrap gap-1 border-0 p-0">
      <legend className="sr-only">Filter by analysis</legend>
      <CommandChip command={null} />
      {commands.map((command) => (
        <CommandChip key={command} command={command} />
      ))}
    </fieldset>
  );
}

function CommandChip({ command }: { command: AppCommand | null }) {
  const commandFilter = useShellStore((state) => state.commandFilter);
  const setCommandFilter = useShellStore((state) => state.setCommandFilter);
  const select = useCallback(() => setCommandFilter(command), [command, setCommandFilter]);

  return (
    <button
      type="button"
      aria-pressed={commandFilter === command}
      onClick={select}
      className="border-line text-ink-muted aria-pressed:bg-accent-subtle aria-pressed:text-ink rounded-full border px-2 py-0.5 text-xs aria-pressed:border-transparent aria-pressed:font-bold"
    >
      {command === null ? "All" : COMMAND_INFO[command].label}
    </button>
  );
}

function RunRow({ run }: { run: RunSummaryResult }) {
  const currentId = useCurrentRunId();
  const compareIds = useShellStore((state) => state.compareIds);
  const toggleCompare = useShellStore((state) => state.toggleCompare);
  const onCompare = useCallback(() => toggleCompare(run.id), [run.id, toggleCompare]);
  const selected = compareIds.includes(run.id);
  const flags = changedFlags(run);
  const headline = run.status === "ok" ? headlineText(run.headline) : "";

  return (
    <div className="relative">
      <Link
        to="/runs/$id/results"
        params={{ id: run.id }}
        aria-current={run.id === currentId ? "page" : undefined}
        className="hover:bg-surface-2 aria-[current=page]:bg-accent-subtle grid grid-cols-[1.125rem_1fr_auto] items-start gap-x-2 gap-y-0.5 rounded-md px-2 py-1.5"
      >
        <span className="mt-0.5" title={statusLabel(run.status)}>
          <StatusIcon status={run.status} size={13} />
        </span>
        <span className="truncate font-bold" title={run.title}>
          {run.title}
        </span>
        <span className="text-ink-muted mt-px text-right text-xs whitespace-nowrap">{headline}</span>
        <span className="text-ink-faint col-span-2 col-start-2 flex items-center gap-1 overflow-hidden pr-6 text-xs">
          <span className="text-accent text-[0.6875rem] font-bold">{run.command}</span>
          {flags.slice(0, CHIPS_SHOWN).map((flag) => (
            <span
              key={flag}
              className="bg-surface-3 text-ink-muted rounded-sm px-1 font-mono text-[0.6875rem] whitespace-nowrap"
            >
              {truncate(flag, CHIP_LENGTH)}
            </span>
          ))}
          {flags.length > CHIPS_SHOWN && <span>+{flags.length - CHIPS_SHOWN}</span>}
        </span>
      </Link>
      {run.status === "ok" && (
        <input
          type="checkbox"
          checked={selected}
          onChange={onCompare}
          title="Select to compare"
          aria-label={`Select ${run.title} to compare`}
          className={cn(
            "absolute right-2 bottom-2 opacity-40 hover:opacity-100 focus-visible:opacity-100",
            selected && "opacity-100",
          )}
        />
      )}
    </div>
  );
}

function CompareBar() {
  const compareIds = useShellStore((state) => state.compareIds);
  const clearCompare = useShellStore((state) => state.clearCompare);
  const [first, second] = compareIds;

  if (first === undefined) {
    return null;
  }

  return (
    <div className="bg-ink text-surface-1 mx-3 mb-3 flex items-center gap-2 rounded-md px-2.5 py-2 text-xs">
      <span className="flex-1">{second === undefined ? "Select one more run to compare" : "2 runs selected"}</span>
      {second !== undefined && (
        <Link
          to="/compare/$a/$b"
          params={{ a: first, b: second }}
          className="border-surface-1/40 rounded-sm border px-2 py-0.5 font-bold"
        >
          Compare
        </Link>
      )}
      <button type="button" onClick={clearCompare} className="px-1 opacity-80 hover:opacity-100">
        Clear
      </button>
    </div>
  );
}
