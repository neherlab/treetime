import type { AppCommand, RunConfigResponseValid, RunSummary } from "@neherlab/app-contracts";

export function duplicateRun(
  runs: readonly RunSummary[],
  command: AppCommand,
  draft: DraftIdentity | undefined,
): DuplicateRun | undefined {
  if (draft?.config_hash === undefined) {
    return undefined;
  }

  const { config_hash: configHash, seed } = draft;

  const run = runs.find(
    (run) =>
      run.status === "ok" &&
      run.command === command &&
      run.config_hash === configHash &&
      (seed === undefined || seed === run.seed),
  );

  return run === undefined ? undefined : { run, repeatsSeed: seed !== undefined };
}

export type DraftIdentity = Pick<RunConfigResponseValid, "config_hash" | "seed">;

export interface DuplicateRun {
  run: RunSummary;
  repeatsSeed: boolean;
}
