import type { RunCheck } from "@neherlab/app-contracts";

export const RUN_BLOCK_TEXT: Readonly<Record<RunBlock, string>> = {
  starting: "Starting the run",
  checking: "Checking the settings",
  blocked: "Resolve the blocking checks first",
  invalid: "Correct the settings marked as invalid first",
};

export function runBlock({ starting, checks, checkPending, valid }: RunReadiness): RunBlock | undefined {
  if (starting) {
    return "starting";
  }

  if (checks === undefined || checkPending) {
    return "checking";
  }

  if (checks.some((check) => check.level === "block")) {
    return "blocked";
  }

  return valid ? undefined : "invalid";
}

export class SingleRun {
  #running = false;

  async run(start: () => Promise<void>): Promise<boolean> {
    if (this.#running) {
      return false;
    }

    this.#running = true;

    try {
      await start();
    } finally {
      this.#running = false;
    }

    return true;
  }
}

export interface RunReadiness {
  starting: boolean;
  checks: readonly RunCheck[] | undefined;
  checkPending: boolean;
  valid: boolean;
}

export type RunBlock = "starting" | "checking" | "blocked" | "invalid";
