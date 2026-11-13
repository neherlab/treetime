import type { JsonValue, UiDraftSource } from "@neherlab/app-contracts";

import { useDraftStore } from "../store/draft";
import { slotHoldingKeys } from "../store/draftRules";
import { releasedUploads } from "./fileSelection";

export function applyDraftInputs(
  assignments: readonly DraftInput[],
  release: (runId: string, paths: readonly string[]) => void,
) {
  const state = useDraftStore.getState();
  const keys = assignments.map(({ key }) => key);

  const draft =
    slotHoldingKeys(state, keys) === "parked" && state.parked !== undefined ? state.parked.draft : state.draft;

  const next = {
    config: { ...draft.config, ...Object.fromEntries(assignments.map(({ key, value }) => [key, value])) },
    sources: { ...draft.sources, ...Object.fromEntries(assignments.map(({ key, source }) => [key, source])) },
  };

  state.assignInputs(next, keys);

  if (draft.upload_run_id !== undefined) {
    release(draft.upload_run_id, releasedUploads(draft, next));
  }
}

export interface DraftInput {
  key: string;
  value: JsonValue;
  source: UiDraftSource;
}
