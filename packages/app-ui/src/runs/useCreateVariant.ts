import { errorMessage, type RunRecord, type SparseConfig } from "@neherlab/app-contracts";
import { runsCreate, runsSettings } from "@neherlab/app-contracts/client";
import { useToastManager } from "@neherlab/design/ui/toast";
import { useNavigate } from "@tanstack/react-router";
import { useCallback, useState } from "react";

import { useApiContext } from "../api/context";
import { commandSettings } from "../settings/catalog";
import { normalizeConfig } from "../settings/config";
import { rerunDraft } from "../settings/rerun";
import { useDraftStore } from "../store/draft";
import { ownDraft, runOptionsDraft } from "../store/draftRules";

export function useCreateVariant(record: RunRecord): () => void {
  const { client } = useApiContext();
  const toasts = useToastManager();

  const create = useCallback(async () => {
    try {
      const { data: settings } = await runsSettings({ client, path: { id: record.id }, throwOnError: true });
      const store = useDraftStore.getState();

      store.load(runOptionsDraft(record, settings, ownDraft(store)));
    } catch (error: unknown) {
      toasts.add({ title: "The settings of the run cannot be loaded", description: errorMessage(error) });
    }
  }, [client, record, toasts]);

  return useCallback(() => void create(), [create]);
}

export function useStartVariant(record: RunRecord, changes: SparseConfig): VariantStart {
  const { client } = useApiContext();
  const navigate = useNavigate();
  const toasts = useToastManager();
  const [busy, setBusy] = useState(false);

  const start = useCallback(async () => {
    setBusy(true);

    try {
      const { data: settings } = await runsSettings({ client, path: { id: record.id }, throwOnError: true });
      const specs = commandSettings(record.command).settings;
      const config = normalizeConfig(specs, { ...rerunDraft(record, settings).config, ...changes });
      const example = record.example === undefined ? {} : { example: record.example };

      const { data: variant } = await runsCreate({
        client,
        body: { command: record.command, config, defer_start: false, source_run_id: record.id, ...example },
        throwOnError: true,
      });

      await navigate({ to: "/runs/$id/results", params: { id: variant.id } });
    } catch (error: unknown) {
      toasts.add({ title: "The variant cannot be started", description: errorMessage(error) });
    } finally {
      setBusy(false);
    }
  }, [changes, client, navigate, record, toasts]);

  return { start: useCallback(() => void start(), [start]), busy };
}

interface VariantStart {
  start: () => void;
  busy: boolean;
}
