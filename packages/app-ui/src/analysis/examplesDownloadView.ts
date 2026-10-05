import type { DatasetCatalog, ExamplesDownloadStatus } from "@neherlab/app-contracts";

export type ExamplesDownloadView =
  | { kind: "hidden" }
  | { kind: "offer" }
  | { kind: "running"; received: number; percent: number | undefined }
  | { kind: "failed"; message: string };

export function examplesDownloadView(
  local: boolean,
  catalog: DatasetCatalog | undefined,
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

  if (!local || catalog === undefined || catalog.datasets.length > 0 || catalog.examples.length > 0) {
    return { kind: "hidden" };
  }

  return download?.state === "failed" ? { kind: "failed", message: download.message } : { kind: "offer" };
}
