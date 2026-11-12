import type { AppCommand, CheckConfigResponse, JobId } from "@neherlab/app-contracts";

const ANALYSIS_CHANGE = "Analysis";

export function staleRun({ shownRunId, runCommand, edited, check }: StaleRunInput): StaleRun {
  if (shownRunId === undefined || runCommand === undefined || check === undefined) {
    return FRESH;
  }

  if (check.command !== runCommand) {
    return { stale: true, changes: [ANALYSIS_CHANGE] };
  }

  if (check.status === "invalid") {
    return edited ? { stale: true, changes: [] } : FRESH;
  }

  const changes = (check.run_differences ?? []).map((difference) => difference.label);

  return changes.length === 0 ? FRESH : { stale: true, changes };
}

export function staleText(changes: readonly string[]): string {
  return changes.length === 0
    ? "Settings changed since this run. Run to update."
    : `Settings changed since this run: ${changesText(changes)}. Run to update.`;
}

export function changesText(changes: readonly string[]): string {
  return changes.map(inSentence).join(", ");
}

export interface StaleRunInput {
  shownRunId: JobId | undefined;
  runCommand: AppCommand | undefined;
  edited: boolean;
  check: CheckConfigResponse | undefined;
}

export type StaleRun = { stale: false } | { stale: true; changes: readonly string[] };

const FRESH: StaleRun = { stale: false };

function inSentence(label: string): string {
  const [first = "", second = ""] = label;

  return second === second.toLowerCase() ? `${first.toLowerCase()}${label.slice(first.length)}` : label;
}
