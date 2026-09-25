import type { AppCommand, RunSummaryResult } from "@neherlab/app-contracts";
import type { DateTime } from "luxon";

import { dayLabel } from "../format";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { COMMAND_INFO } from "../settings/commands";
import { settingFlags } from "../settings/titles";

interface RunGroup {
  label: string;
  runs: RunSummaryResult[];
}

export function listedRuns(
  runs: readonly RunSummaryResult[],
  filter: string,
  command: AppCommand | null,
): RunSummaryResult[] {
  const words = filter
    .toLowerCase()
    .split(/\s+/u)
    .filter((word) => word !== "");

  return runs.filter((run) => {
    if (run.status === "created" || (command !== null && run.command !== command)) {
      return false;
    }

    const haystack = [run.title, run.command, COMMAND_INFO[run.command].label, ...changedFlags(run)]
      .join(" ")
      .toLowerCase();

    return words.every((word) => haystack.includes(word));
  });
}

export function groupRuns(runs: readonly RunSummaryResult[], now: DateTime): RunGroup[] {
  const groups: RunGroup[] = [];
  const pinned = runs.filter((run) => run.pinned);

  if (pinned.length > 0) {
    groups.push({ label: "Pinned", runs: pinned });
  }

  for (const run of runs) {
    if (run.pinned) {
      continue;
    }

    const label = dayLabel(run.created_at, now);
    const group = groups.find((candidate) => candidate.label === label && candidate.label !== "Pinned");

    if (group === undefined) {
      groups.push({ label, runs: [run] });
    } else {
      group.runs.push(run);
    }
  }

  return groups;
}

export function changedFlags(run: RunSummaryResult): string[] {
  return settingFlags(COMMAND_SETTINGS[run.command].specs, run.changed_settings);
}
