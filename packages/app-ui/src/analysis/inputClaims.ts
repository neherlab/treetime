import { useDraftStore } from "../store/draft";

export function claimInputs(keys: readonly string[]): InputClaim {
  const selection = Symbol("input claim");
  const { epoch, claimSelections } = useDraftStore.getState();
  claimSelections(keys, selection);

  return (runId) => {
    const { draft, epoch: current, selections } = useDraftStore.getState();

    return (
      current === epoch &&
      keys.every((key) => selections[key] === selection) &&
      (runId === undefined || draft.upload_run_id === runId)
    );
  };
}

export type InputClaim = (runId?: string) => boolean;
