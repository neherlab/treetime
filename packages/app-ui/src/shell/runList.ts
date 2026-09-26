import type { AppCommand, RunSummary } from "@neherlab/app-contracts";
import type { DateTime } from "luxon";

import { dayLabel } from "../format";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { COMMAND_INFO } from "../settings/commands";
import { settingFlags } from "../settings/titles";
import { wordMatcher } from "../text";

interface RunGroup {
  label: string;
  runs: RunSummary[];
}

export function listedRuns(runs: readonly RunSummary[], filter: string, command: AppCommand | null): RunSummary[] {
  const matches = wordMatcher(filter);

  return runs.filter(
    (run) =>
      run.status !== "created" &&
      (command === null || run.command === command) &&
      matches([run.title, run.command, COMMAND_INFO[run.command].label, ...changedFlags(run)].join(" ")),
  );
}

export function groupRuns(runs: readonly RunSummary[], now: DateTime): RunGroup[] {
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

export function changedFlags(run: RunSummary): string[] {
  return settingFlags(COMMAND_SETTINGS[run.command].specs, run.changed_settings);
}
