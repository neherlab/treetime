import type { ExamplesDownloadStatus } from "@neherlab/app-contracts";
import { datasets, examplesDownload, examplesDownloadStart } from "@neherlab/app-contracts/client";
import { useCallback } from "react";
import Download from "~icons/lucide/download";

import { useApi, useApiMutation } from "../api/hooks";
import { formatBytes } from "../format";
import { useHost } from "../host-context";
import { newerDownloadStatus, useExamplesDownloadStore } from "../store/examplesDownload";
import { Alert, AlertDescription } from "../ui/alert";
import { Button } from "../ui/button";
import { Progress } from "../ui/progress";
import { Spinner } from "../ui/spinner";
import { examplesDownloadView } from "./examplesDownloadView";

export function ExamplesDownload() {
  const host = useHost();
  const { data: catalog } = useApi((context) => datasets(context), { staleTime: Infinity });
  const status = useExamplesDownloadStatus(host !== null);
  const view = examplesDownloadView(host !== null, catalog, status);
  const { mutate: start, isPending, error } = useStartExamplesDownload();
  const onStart = useCallback(() => start(undefined), [start]);

  if (view.kind === "hidden") {
    return null;
  }

  if (view.kind === "running") {
    return (
      <div className="grid gap-1.5 px-3.5 py-3 text-sm">
        <span className="text-muted-foreground">Downloading the example datasets: {formatBytes(view.received)}</span>
        <Progress value={view.percent ?? null} aria-label="Download of the example datasets" />
      </div>
    );
  }

  return (
    <div className="grid gap-2 px-3.5 py-3 text-sm">
      <p className="text-muted-foreground">
        The examples folder is empty. Download the example datasets and configs of this TreeTime release.
      </p>
      {view.kind === "failed" && (
        <Alert variant="destructive">
          <AlertDescription>The download failed: {view.message}</AlertDescription>
        </Alert>
      )}
      {error !== null && (
        <Alert variant="destructive">
          <AlertDescription>The download cannot start: {error.message}</AlertDescription>
        </Alert>
      )}
      <div>
        <Button type="button" variant="outline" size="sm" disabled={isPending} onClick={onStart}>
          {isPending ? <Spinner /> : <Download aria-hidden />}
          {view.kind === "failed" ? "Try again" : "Download examples"}
        </Button>
      </div>
    </div>
  );
}

export function useExamplesDownloadStatus(local: boolean): ExamplesDownloadStatus | undefined {
  const { data } = useApi((context) => examplesDownload(context), { enabled: local });
  const fromEvents = useExamplesDownloadStore((state) => state.status);

  return newerDownloadStatus(fromEvents, data);
}

export function useStartExamplesDownload() {
  return useApiMutation((context) => examplesDownloadStart(context), {
    seed: () => (context) => examplesDownload(context),
  });
}
