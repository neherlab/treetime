import { type SparseConfig, errorMessage, type AppCommand } from "@neherlab/app-contracts";
import { runsCreate, runsStart } from "@neherlab/app-contracts/client";
import { useNavigate } from "@tanstack/react-router";
import { useCallback } from "react";

import { useApiContext } from "../api/context";
import { useDraftStore } from "../store/draft";
import { useToastManager } from "../ui/toast";
import { pendingUploadRun } from "./pendingUpload";

export function useStartRun(command: AppCommand) {
  const { client } = useApiContext();
  const navigate = useNavigate();
  const toasts = useToastManager();

  return useCallback(
    async (config: SparseConfig) => {
      const { draft, clearRuns } = useDraftStore.getState();

      try {
        const upload = await pendingUploadRun(client, draft.upload_run_id);

        const { data: record } =
          upload === undefined
            ? await runsCreate({ client, body: { command, config, defer_start: false }, throwOnError: true })
            : await runsStart({ client, path: { id: upload.id }, body: { command, config }, throwOnError: true });

        clearRuns();
        await navigate({ to: "/runs/$id/results", params: { id: record.id } });
      } catch (error: unknown) {
        toasts.add({
          title: "The run cannot be started",
          description: errorMessage(error),
        });
      }
    },
    [client, command, navigate, toasts],
  );
}
