import type { RunRecord } from "@neherlab/app-contracts";
import { useNavigate } from "@tanstack/react-router";
import { useCallback } from "react";

import { rerunDraft } from "../settings/rerun";
import { useDraftStore } from "../store/draft";

export function useRerun(record: RunRecord): () => void {
  const navigate = useNavigate();

  return useCallback(() => {
    const draft = rerunDraft(record);

    useDraftStore.getState().load({
      command: record.command,
      config: draft.config,
      sources: Object.fromEntries(
        Object.entries(draft.inputLabels).map(([key, label]) => [key, { label, origin: "run" as const, size: null }]),
      ),
      from_run_id: record.id,
      upload_run_id: null,
    });
    void navigate({ to: "/new" });
  }, [navigate, record]);
}
