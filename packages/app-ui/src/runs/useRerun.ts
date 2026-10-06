import type { AppCommand, RunRecord } from "@neherlab/app-contracts";
import { useNavigate } from "@tanstack/react-router";
import { useCallback } from "react";

import { rerunDraft } from "../settings/rerun";
import { useDraftStore } from "../store/draft";

export function useRerun(record: RunRecord, command: AppCommand = record.command): () => void {
  const navigate = useNavigate();

  return useCallback(() => {
    const draft = rerunDraft(record, command);

    const { draft: current, load } = useDraftStore.getState();
    const { upload_run_id: _uploadRun, ...kept } = current;

    load({
      ...kept,
      command,
      config: draft.config,
      sources: Object.fromEntries(
        Object.entries(draft.inputLabels).map(([key, label]) => [key, { label, origin: "run" as const }]),
      ),
      from_run_id: record.id,
    });
    void navigate({ to: "/new" });
  }, [command, navigate, record]);
}
