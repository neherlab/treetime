import { errorMessage } from "@neherlab/app-contracts";
import { Alert, AlertDescription, AlertTitle } from "@neherlab/design/ui/alert";
import { Button } from "@neherlab/design/ui/button";
import { toast } from "@neherlab/design/ui/toast";
import { useQueryClient } from "@tanstack/react-query";
import { useCallback, useEffect } from "react";

import { resetApiQueries } from "../api/hooks";
import type { Host } from "../host";
import { useHost } from "../host-context";
import { type BackendStopNotice, stopAction, stopEffect, useAppNoticesStore } from "./noticeStore";

export function AppNotices() {
  const host = useHost();
  const stop = useAppNoticesStore((state) => state.stop);
  const stops = useAppNoticesStore((state) => state.stops);
  const chunkFailed = useAppNoticesStore((state) => state.chunkFailed);
  const queryClient = useQueryClient();

  useEffect(() => {
    if (stops > 0) {
      void resetApiQueries(queryClient);
    }
  }, [queryClient, stops]);

  if (stop === undefined && !chunkFailed) {
    return null;
  }

  return (
    <div className="grid gap-2 px-6 pt-4">
      {stop !== undefined && <StopBanner stop={stop} host={host} />}
      {chunkFailed && (
        <Alert variant="warning">
          <AlertTitle>A part of TreeTime cannot be loaded</AlertTitle>
          <AlertDescription>
            A newer version of TreeTime was published, or the network failed. Reloading the page loads TreeTime again
            and loses the runs of this tab; download run packages first to keep runs.
          </AlertDescription>
          <NoticeActions>
            <Button type="button" variant="outline" size="sm" onClick={reloadPage}>
              Reload the page
            </Button>
          </NoticeActions>
        </Alert>
      )}
    </div>
  );
}

function StopBanner({ stop, host }: { stop: BackendStopNotice; host: Host | null }) {
  const dismissStop = useAppNoticesStore((state) => state.dismissStop);
  const action = stopAction(stop.restarts, host !== null);

  const restart = useCallback(() => {
    if (host !== null) {
      void restartBackEnd(host, dismissStop);
    }
  }, [dismissStop, host]);

  return (
    <Alert variant="destructive">
      <AlertTitle>The back end stopped</AlertTitle>
      <AlertDescription>
        <p>{stop.message}</p>
        <p>{stopEffect(host !== null)}</p>
      </AlertDescription>
      <NoticeActions>
        {action === "dismiss" && (
          <Button type="button" variant="outline" size="sm" onClick={dismissStop}>
            Dismiss
          </Button>
        )}
        {action === "reload-page" && (
          <Button type="button" variant="outline" size="sm" onClick={reloadPage}>
            Reload the page
          </Button>
        )}
        {action === "restart-back-end" && (
          <Button type="button" variant="outline" size="sm" onClick={restart}>
            Restart the back end
          </Button>
        )}
      </NoticeActions>
    </Alert>
  );
}

async function restartBackEnd(host: Host, dismissStop: () => void): Promise<void> {
  try {
    await host.restartBackend();
    dismissStop();
  } catch (error: unknown) {
    toast.add({ title: "The back end cannot be restarted", description: errorMessage(error) });
  }
}

function NoticeActions({ children }: { children: React.ReactNode }) {
  return <div className="mt-2 flex gap-2">{children}</div>;
}

function reloadPage(): void {
  window.location.reload();
}
