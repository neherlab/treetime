import { errorMessage, type RunRecord } from "@neherlab/app-contracts";
import { runsSettings } from "@neherlab/app-contracts/client";
import { useToastManager } from "@neherlab/design/ui/toast";
import { useNavigate } from "@tanstack/react-router";
import { useCallback } from "react";

import { useApiContext } from "../api/context";
import { rerunDraft } from "../settings/rerun";
import { useDraftStore } from "../store/draft";

export function useRerun(record: RunRecord): () => void {
  const { client } = useApiContext();
  const navigate = useNavigate();
  const toasts = useToastManager();

  const rerun = useCallback(async () => {
    try {
      const { data: settings } = await runsSettings({ client, path: { id: record.id }, throwOnError: true });
      const draft = rerunDraft(record, settings);

      const { draft: current, load } = useDraftStore.getState();
      const { upload_run_id: _uploadRun, ...kept } = current;

      load({
        ...kept,
        command: record.command,
        config: draft.config,
        sources: Object.fromEntries(
          Object.entries(draft.inputLabels).map(([key, label]) => [
            key,
            { label, origin: { kind: "run" as const }, run_id: record.id },
          ]),
        ),
        from_run_id: record.id,
      });
      await navigate({ to: "/new" });
    } catch (error: unknown) {
      toasts.add({ title: "The settings of the run cannot be loaded", description: errorMessage(error) });
    }
  }, [client, navigate, record, toasts]);

  return useCallback(() => void rerun(), [rerun]);
}
