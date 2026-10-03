import type { AppCommand, UiDraftSource } from "@neherlab/app-contracts";
import { create } from "zustand";

import type { JsonObject } from "../settings/json";
import { freshDraft, type Draft } from "./draftSchema";

interface DraftState extends Draft {
  epoch: number;
}

interface DraftActions {
  load: (draft: Partial<Draft>) => void;
  update: (draft: Partial<Omit<Draft, "command" | "config">>) => void;
  setConfig: (config: JsonObject) => void;
  setSource: (key: string, source: UiDraftSource | null) => void;
  reset: (command: AppCommand) => void;
}

export const useDraftStore = create<DraftState & DraftActions>()((set) => ({
  ...freshDraft("timetree"),
  epoch: 0,
  load: (draft) => {
    set((state) => ({ ...draft, epoch: state.epoch + 1 }));
  },
  update: (draft) => {
    set(draft);
  },
  setConfig: (config) => {
    set({ config });
  },
  setSource: (key, source) => {
    set((state) => {
      const sources = Object.fromEntries(Object.entries(state.sources).filter(([name]) => name !== key));

      return { sources: source === null ? sources : { ...sources, [key]: source } };
    });
  },
  reset: (command) => {
    set((state) => ({ ...freshDraft(command), epoch: state.epoch + 1 }));
  },
}));
