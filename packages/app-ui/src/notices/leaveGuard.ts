import type { UiDraft } from "@neherlab/app-contracts";
import { create } from "zustand";

export interface TabWork {
  runs: number;
  imports: number;
  formUploads: number;
}

export function leaveLosesWork({ runs, imports, formUploads }: TabWork): boolean {
  return runs > 0 || imports > 0 || formUploads > 0;
}

export function formUploads(draft: UiDraft): number {
  return Object.values(draft.sources).filter((source) => source.origin.kind === "upload").length;
}

interface TabWorkState {
  imports: number;
  holdsWork: boolean;
}

interface TabWorkActions {
  startImport: () => () => void;
  setHoldsWork: (holdsWork: boolean) => void;
}

export const useTabWorkStore = create<TabWorkState & TabWorkActions>()((set) => ({
  imports: 0,
  holdsWork: false,
  startImport: () => {
    set(({ imports }) => ({ imports: imports + 1 }));
    let ended = false;

    return () => {
      if (!ended) {
        ended = true;
        set(({ imports }) => ({ imports: imports - 1 }));
      }
    };
  },
  setHoldsWork: (holdsWork) => set({ holdsWork }),
}));
