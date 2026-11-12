import { useDraftStore } from "../store/draft";
import { claimedSlot } from "../store/draftRules";

export function claimInputs(keys: readonly string[]): InputClaim {
  const selection = Symbol("input claim");
  useDraftStore.getState().claimSelections(keys, selection);

  const isCurrent = (runId?: string) => {
    const state = useDraftStore.getState();
    const slot = claimedSlot(state, selection) === "parked" ? state.parked : undefined;
    const { draft, selections } = slot ?? state;

    return keys.every((key) => selections[key] === selection) && (runId === undefined || draft.upload_run_id === runId);
  };

  return Object.assign(isCurrent, { selection });
}

export type InputClaim = ((runId?: string) => boolean) & { selection: symbol };
