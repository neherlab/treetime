<<<<<<< HEAD
<<<<<<< HEAD
export interface CommandResult {
  status: "ok" | "error";
  error?: string;
=======
export class CancelledError extends Error {
  constructor() {
    super("Operation cancelled");
    this.name = "CancelledError";
  }
>>>>>>> ac719231 (feat(web): wire AbortController through bridge contract and UI)
}

=======
>>>>>>> c2b9da5e (feat: add per-command result types and hooks across TypeScript layer)
export interface ProgressEvent {
  stage: string;
  fraction: number;
  message: string;
}

export interface ErrorResponse {
  code: string;
  message: string;
}

export interface AncestralResult {
  model_name: string;
}

export interface ClockResult {
  clock_model: Record<string, unknown>;
  regression_results: Record<string, unknown>[];
}

export interface TimetreeResult {
  clock_model: Record<string, unknown>;
}

export interface MugrationResult {
  log_lh: number;
}

export interface OptimizeResult {}

export interface PruneResult {}

export interface HomoplasyResult {}
