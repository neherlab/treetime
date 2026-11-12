import type { JobId, RunRecord, SparseConfig } from "@neherlab/app-contracts";
import { runsGet, runsSettings } from "@neherlab/app-contracts/client";
import { useParams } from "@tanstack/react-router";
import { useEffect } from "react";

import { useApi } from "../api/hooks";
import { LoadingState } from "../components/PageShell";
import { useDraftStore } from "../store/draft";
import { ownDraft, runOptionsDraft } from "../store/draftRules";
import { DraftForm } from "./DraftForm";

export function WorkspaceRoute() {
  const { id: shownRunId } = useParams({ strict: false });
  const pending = useShownRunOptions(shownRunId);
  const command = useDraftStore((state) => state.draft.command);
  const epoch = useDraftStore((state) => state.epoch);

  if (pending) {
    return <LoadingState text="Loading the options of the run" />;
  }

  return <DraftForm key={`${command}:${epoch}:${shownRunId ?? "new"}`} command={command} shownRunId={shownRunId} />;
}

function useShownRunOptions(shownRunId: JobId | undefined): boolean {
  const enabled = shownRunId !== undefined;
  const id = shownRunId ?? "";

  const { data: record, error: recordError } = useApi((context) => runsGet({ ...context, path: { id } }), {
    enabled,
  });

  const { data: settings, error: settingsError } = useApi((context) => runsSettings({ ...context, path: { id } }), {
    enabled,
    staleTime: Infinity,
  });

  const lastShown = useDraftStore((state) => state.lastShown);

  useEffect(() => {
    showRunOptions(shownRunId, record, settings);
  }, [record, settings, shownRunId]);

  return enabled && lastShown !== shownRunId && recordError === null && settingsError === null;
}

function showRunOptions(
  shownRunId: JobId | undefined,
  record: RunRecord | undefined,
  settings: SparseConfig | undefined,
): void {
  const store = useDraftStore.getState();

  if (shownRunId === undefined) {
    store.leaveRun();

    return;
  }

  if (record?.id !== shownRunId || settings === undefined) {
    return;
  }

  const options = runOptionsDraft(record, settings, ownDraft(store));

  if (store.lastShown === shownRunId) {
    store.refreshView(shownRunId, options);
  } else {
    store.viewRun(shownRunId, options);
  }
}
