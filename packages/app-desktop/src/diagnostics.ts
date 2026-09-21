import { mkdirSync } from "node:fs";
import { resolve } from "node:path";
import { setFlagsFromString } from "node:v8";

export function initDiagnostics(title: string): void {
  process.title = title;

  Error.stackTraceLimit = 50;
  process.setSourceMapsEnabled(true);

  process.on("warning", (warning) => {
    console.warn(warning.stack ?? `${warning.name}: ${warning.message}`);
  });

  const diagnosticDir = resolve(process.env["TREETIME_DIAGNOSTIC_DIR"] ?? resolve(process.cwd(), "tmp/diagnostics"));
  mkdirSync(diagnosticDir, { recursive: true });

  if (process.report) {
    process.report.directory = diagnosticDir;
    process.report.reportOnFatalError = true;
    process.report.reportOnSignal = true;
    process.report.signal = "SIGUSR2";
  }

  setFlagsFromString("--heapsnapshot-near-heap-limit=1");
}
