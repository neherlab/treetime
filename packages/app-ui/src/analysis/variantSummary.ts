import type { CheckConfigResponse, JobId, RunSummary, SettingSpec } from "@neherlab/app-contracts";

import { changesText, staleRun } from "./staleRun";

export function variantSummary({
  fromRunId,
  shownRunId,
  viewing,
  source,
  check,
  specs,
}: VariantInput): string | undefined {
  if (viewing || fromRunId === undefined) {
    return undefined;
  }

  const onSource = fromRunId === shownRunId;
  const name = onSource ? "this run" : `"${source?.title ?? fromRunId}"`;
  const comparedWithSource = shownRunId === undefined || onSource;

  if (!comparedWithSource || source === undefined || check === undefined || check.status === "invalid") {
    return `Variant of ${name}`;
  }

  const stale = staleRun({ shownRunId: fromRunId, runCommand: source.command, edited: true, check, specs });

  if (!stale.stale) {
    return `Variant of ${name}. No settings changed yet.`;
  }

  return onSource ? `Variant of ${name}` : `Variant of ${name}. Changed: ${changesText(stale.changes)}.`;
}

export interface VariantInput {
  fromRunId: JobId | undefined;
  shownRunId: JobId | undefined;
  viewing: boolean;
  source: Pick<RunSummary, "title" | "command"> | undefined;
  check: CheckConfigResponse | undefined;
  specs: readonly SettingSpec[];
}
