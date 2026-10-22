import { useDraftStore } from "../store/draft";

export function claimInputs(keys: readonly string[]): InputClaim {
  const selection = Symbol("input claim");
  useDraftStore.getState().claimSelections(keys, selection);

  return (runId) => {
    const { draft, selections } = useDraftStore.getState();

    return keys.every((key) => selections[key] === selection) && (runId === undefined || draft.upload_run_id === runId);
  };
}

export type InputClaim = (runId?: string) => boolean;
