import { errorMessage } from "@neherlab/app-contracts";
import type { AppCommand } from "@neherlab/app-contracts";
import { useQueryClient } from "@tanstack/react-query";
import { useNavigate } from "@tanstack/react-router";
import { useCallback } from "react";

import { useBridge } from "../BridgeContext";
import { RUNS_KEY } from "../queries";
import type { JsonObject } from "../settings/json";
import { useDraftStore } from "../store/draft";
import { Toast } from "../ui";
import { pendingUploadRun } from "./pendingUpload";

export function useStartRun(command: AppCommand) {
  const bridge = useBridge();
  const navigate = useNavigate();
  const queryClient = useQueryClient();
  const toasts = Toast.useToastManager();

  return useCallback(
    async (config: JsonObject, fallbackTitle: string) => {
      const draft = useDraftStore.getState();

      const title = draft.title.trim() === "" ? fallbackTitle : draft.title.trim();

      try {
        const upload = await pendingUploadRun(bridge, draft.uploadRunId);
        let id: string;

        if (upload !== undefined && upload.command === command) {
          await bridge.updateRun(upload.id, { title });
          id = (await bridge.startRun(upload.id, { config })).id;
        } else {
          id = (await bridge.createRun({ command, config, title, defer_start: false })).id;
        }

        draft.update({ title: "", fromRunId: null, uploadRunId: null });
        await queryClient.invalidateQueries({ queryKey: RUNS_KEY });
        await navigate({ to: "/runs/$id/results", params: { id } });
      } catch (error: unknown) {
        toasts.add({
          title: "The run cannot be started",
          description: errorMessage(error),
        });
      }
    },
    [bridge, command, navigate, queryClient, toasts],
  );
}
