import type { JobId, RunRecord, SparseConfig } from "@neherlab/app-contracts";
import { useParams } from "@tanstack/react-router";
import { useEffect } from "react";

import { useDraftStore } from "../store/draft";
import { ownDraft, runOptionsDraft } from "../store/draftRules";
import { DraftForm } from "./DraftForm";
import { useShownRunSettings } from "./shownRun";

export function WorkspaceRoute() {
  const { id: shownRunId } = useParams({ strict: false });
  const pending = useShownRunOptions(shownRunId);
  const command = useDraftStore((state) => state.draft.command);

  return <DraftForm command={command} shownRunId={shownRunId} pending={pending} />;
}

function useShownRunOptions(shownRunId: JobId | undefined): boolean {
  const { record, settings } = useShownRunSettings(shownRunId);
  const lastShown = useDraftStore((state) => state.lastShown);

  useEffect(() => {
    showRunOptions(shownRunId, record.data, settings.data);
  }, [record.data, settings.data, shownRunId]);

  return shownRunId !== undefined && lastShown !== shownRunId && record.error === null && settings.error === null;
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

  if (options === undefined) {
    return;
  }

  if (store.lastShown === shownRunId) {
    store.refreshView(shownRunId, options);
  } else {
    store.viewRun(shownRunId, options);
  }
}
