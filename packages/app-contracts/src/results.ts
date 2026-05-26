export interface CommandResult {
  status: "ok" | "error";
  error?: string;
}

export interface ProgressEvent {
  stage: string;
  fraction: number;
  message: string;
}
