import type { ExamplesDownloadStatus } from "@neherlab/app-contracts";

import type { ExampleIndex } from "./exampleData";

export type ExamplesDownloadView =
  | { kind: "hidden" }
  | { kind: "offer" }
  | { kind: "running"; received: number; percent: number | undefined }
  | { kind: "failed"; message: string };

export function examplesDownloadView(
  examples: ExampleIndex | undefined,
  status: ExamplesDownloadStatus | undefined,
): ExamplesDownloadView {
  const download = status?.download;

  if (download?.state === "running") {
    const total = download.total ?? 0;

    return {
      kind: "running",
      received: download.received,
      percent: total > 0 ? Math.min(100, Math.round((download.received / total) * 100)) : undefined,
    };
  }

  if (examples?.origin !== "local" || examples.index !== undefined) {
    return { kind: "hidden" };
  }

  return download?.state === "failed" ? { kind: "failed", message: download.message } : { kind: "offer" };
}
