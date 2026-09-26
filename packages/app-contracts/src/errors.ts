export class CancelledError extends Error {
  constructor() {
    super("Operation cancelled");
    this.name = "CancelledError";
  }
}

export class CommandError extends Error {
  readonly jobId: string;
  readonly causes: string[];

  constructor(jobId: string, message: string, causes: string[]) {
    super(message);
    this.name = "CommandError";
    this.jobId = jobId;
    this.causes = causes;
  }
}

export class RunEndedError extends Error {
  readonly runId: string;

  constructor(runId: string) {
    super(`the event stream of run ${runId} ended without a terminal event`);
    this.name = "RunEndedError";
    this.runId = runId;
  }
}

export function errorMessage(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
