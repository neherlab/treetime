import type { AppCommand, JobId, RunSummary } from "@neherlab/app-contracts";
import { DateTime } from "luxon";

import { fuzzyFilter } from "../fuzzy";
import { commandSettings } from "../settings/catalog";

export interface RunGroup {
  id: JobId;
  entries: RunEntry[];
}

export interface RunEntry {
  run: RunSummary;
  depth: number;
}

export function listedRuns(
  runs: readonly RunSummary[],
  filter: string,
  command: AppCommand | null,
  hidden: ReadonlySet<JobId> = new Set(),
): RunSummary[] {
  const candidates = runs.filter(
    (run) => run.status !== "created" && !hidden.has(run.id) && (command === null || run.command === command),
  );

  return fuzzyFilter(candidates, filter, (run) =>
    [run.title, run.command, commandSettings(run.command).title].join(" "),
  );
}

export function groupRuns(runs: readonly RunSummary[]): RunGroup[] {
  const ids = new Set(runs.map((run) => run.id));
  const byGroup = Object.groupBy(runs, (run) => run.group_id);

  return Object.entries(byGroup)
    .flatMap(([id, members]) => (members === undefined ? [] : [{ id, members }]))
    .toSorted((left, right) => compareGroups(right.members, left.members))
    .map(({ id, members }) => ({ id, entries: derivationOrder(members, ids) }));
}

export function newestRun(group: RunGroup): RunSummary | undefined {
  return group.entries.map((entry) => entry.run).toSorted((left, right) => created(right) - created(left))[0];
}

export function groupedRunOrder(groups: readonly RunGroup[]): RunSummary[] {
  return groups.flatMap((group) => group.entries.map((entry) => entry.run));
}

function derivationOrder(members: readonly RunSummary[], listed: ReadonlySet<JobId>): RunEntry[] {
  const sorted = members.toSorted((left, right) => created(left) - created(right));
  const earlier = new Map(sorted.map((run, index) => [run.id, index]));

  const parentOf = (run: RunSummary, index: number) => {
    const source = run.source_run_id;

    return source !== undefined && listed.has(source) && (earlier.get(source) ?? index) < index ? source : "";
  };

  const children = Object.groupBy(sorted, parentOf);

  const visit = (run: RunSummary, depth: number): RunEntry[] => [
    { run, depth },
    ...(children[run.id] ?? []).flatMap((child) => visit(child, depth + 1)),
  ];

  return (children[""] ?? []).flatMap((root) => visit(root, 0));
}

function compareGroups(left: readonly RunSummary[], right: readonly RunSummary[]): number {
  return Number(isPinned(left)) - Number(isPinned(right)) || latest(left) - latest(right);
}

function isPinned(members: readonly RunSummary[]): boolean {
  return members.some((run) => run.pinned);
}

function latest(members: readonly RunSummary[]): number {
  return Math.max(...members.map(created));
}

function created(run: RunSummary): number {
  return DateTime.fromISO(run.created_at).toMillis();
}
