import { errorMessage, type RunRecord } from "@neherlab/app-contracts";
import { runsSettings } from "@neherlab/app-contracts/client";
import { useToastManager } from "@neherlab/design/ui/toast";
import { useNavigate } from "@tanstack/react-router";
import { useCallback } from "react";

import { useApiContext } from "../api/context";
import { useDraftStore } from "../store/draft";
import { ownDraft, runOptionsDraft } from "../store/draftRules";

export function useRerun(record: RunRecord): () => void {
  const { client } = useApiContext();
  const navigate = useNavigate();
  const toasts = useToastManager();

  const rerun = useCallback(async () => {
    try {
      const { data: settings } = await runsSettings({ client, path: { id: record.id }, throwOnError: true });
      const store = useDraftStore.getState();
      const options = runOptionsDraft(record, settings, ownDraft(store));

      if (options === undefined) {
        throw new Error("the run has not resolved its settings yet");
      }

      store.load(options);
      await navigate({ to: "/new" });

      const [named] = record.error?.settings ?? [];

      if (named !== undefined) {
        useDraftStore.getState().focusSetting(named);
      }
    } catch (error: unknown) {
      toasts.add({ title: "The settings of the run cannot be loaded", description: errorMessage(error) });
    }
  }, [client, navigate, record, toasts]);

  return useCallback(() => void rerun(), [rerun]);
}
