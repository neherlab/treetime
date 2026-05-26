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

export interface ProgressEvent {
  stage: string;
  fraction: number;
  message: string;
}
