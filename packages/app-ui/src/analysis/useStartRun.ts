import { errorMessage } from "@neherlab/app-contracts";
import type { AppCommand } from "@neherlab/app-contracts";
import { runsCreate, runsStart } from "@neherlab/app-contracts/client";
import { useNavigate } from "@tanstack/react-router";
import { useCallback } from "react";

import { useApiContext } from "../api/context";
import type { JsonObject } from "../settings/json";
import { useDraftStore } from "../store/draft";
import { Toast } from "../ui";
import { pendingUploadRun } from "./pendingUpload";

export function useStartRun(command: AppCommand) {
  const { client } = useApiContext();
  const navigate = useNavigate();
  const toasts = Toast.useToastManager();

  return useCallback(
    async (config: JsonObject) => {
      const draft = useDraftStore.getState();

      try {
        const upload = await pendingUploadRun(client, draft.uploadRunId);
        let id: string;

        if (upload !== undefined && upload.command === command) {
          id = (await runsStart({ client, path: { id: upload.id }, body: { config }, throwOnError: true })).data.id;
        } else {
          id = (await runsCreate({ client, body: { command, config, defer_start: false }, throwOnError: true })).data
            .id;
        }

        draft.update({ fromRunId: null, uploadRunId: null });
        await navigate({ to: "/runs/$id/results", params: { id } });
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
