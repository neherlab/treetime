import type { JobId, RunSummary } from "@neherlab/app-contracts";

export interface ComparePair {
  a: JobId;
  b: JobId;
}

export function comparePair(
  selectedIds: readonly JobId[],
  runs: ReadonlyArray<Pick<RunSummary, "id" | "status">>,
): ComparePair | undefined {
  const [a, b, ...others] = selectedIds;

  if (a === undefined || b === undefined || others.length > 0) {
    return undefined;
  }

  const finished = (id: JobId) => runs.some((run) => run.id === id && run.status === "ok");

  return finished(a) && finished(b) ? { a, b } : undefined;
}
