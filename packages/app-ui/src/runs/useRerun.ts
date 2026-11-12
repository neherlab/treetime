import { errorMessage, type RunRecord } from "@neherlab/app-contracts";
import { runsSettings } from "@neherlab/app-contracts/client";
import { useToastManager } from "@neherlab/design/ui/toast";
import { useNavigate } from "@tanstack/react-router";
import { useCallback } from "react";

import { useApiContext } from "../api/context";
import { useDraftStore } from "../store/draft";
import { runOptionsDraft } from "../store/draftRules";

export function useRerun(record: RunRecord): () => void {
  const { client } = useApiContext();
  const navigate = useNavigate();
  const toasts = useToastManager();

  const rerun = useCallback(async () => {
    try {
      const { data: settings } = await runsSettings({ client, path: { id: record.id }, throwOnError: true });
      const { draft, viewing, load } = useDraftStore.getState();

      load(runOptionsDraft(record, settings, viewing?.parked ?? draft));
      await navigate({ to: "/new" });
    } catch (error: unknown) {
      toasts.add({ title: "The settings of the run cannot be loaded", description: errorMessage(error) });
    }
  }, [client, navigate, record, toasts]);

  return useCallback(() => void rerun(), [rerun]);
}
