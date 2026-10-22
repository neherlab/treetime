import type { JsonValue, UiDraftSource } from "@neherlab/app-contracts";

import { useDraftStore } from "../store/draft";
import { releasedUploads } from "./fileSelection";

export function applyDraftInputs(
  assignments: readonly DraftInput[],
  release: (runId: string, paths: readonly string[]) => void,
) {
  const { draft, assignInputs } = useDraftStore.getState();

  const next = {
    config: { ...draft.config, ...Object.fromEntries(assignments.map(({ key, value }) => [key, value])) },
    sources: { ...draft.sources, ...Object.fromEntries(assignments.map(({ key, source }) => [key, source])) },
  };

  assignInputs(
    next,
    assignments.map(({ key }) => key),
  );

  if (draft.upload_run_id !== undefined) {
    release(draft.upload_run_id, releasedUploads(draft, next));
  }
}

export interface DraftInput {
  key: string;
  value: JsonValue;
  source: UiDraftSource;
}
