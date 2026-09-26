import { mkdirSync } from "node:fs";
import { setHeapSnapshotNearHeapLimit } from "node:v8";

export const DIAGNOSTIC_DIR_ENV = "TREETIME_DIAGNOSTIC_DIR";

export function initDiagnostics(title: string, diagnosticDir: string): void {
  process.title = title;

  Error.stackTraceLimit = 50;
  process.setSourceMapsEnabled(true);

  process.on("warning", (warning) => {
    console.warn(warning.stack ?? `${warning.name}: ${warning.message}`);
  });

  mkdirSync(diagnosticDir, { recursive: true });

  if (process.report !== undefined) {
    process.report.directory = diagnosticDir;
    process.report.reportOnFatalError = true;
    process.report.reportOnSignal = true;
    process.report.signal = "SIGUSR2";
  }

  setHeapSnapshotNearHeapLimit(1);
}
