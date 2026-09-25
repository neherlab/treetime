import type { AppCommand } from "@neherlab/app-contracts";
import { create } from "zustand";
import { createJSONStorage, persist } from "zustand/middleware";

import type { JsonObject } from "../settings/json";
import { freshDraft, restoredDraft, zDraft, zStoredValues, type Draft, type InputSource } from "./draftSchema";

const DRAFT_STORAGE_KEY = "treetime-draft";

const DRAFT_VERSION = 1;

interface DraftState extends Draft {
  epoch: number;
}

interface DraftActions {
  load: (draft: Partial<Draft>) => void;
  update: (draft: Partial<Omit<Draft, "command" | "config">>) => void;
  setConfig: (config: JsonObject) => void;
  setSource: (key: string, source: InputSource | null) => void;
  reset: (command: AppCommand) => void;
}

export const useDraftStore = create<DraftState & DraftActions>()(
  persist(
    (set) => ({
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
    }),
    {
      name: DRAFT_STORAGE_KEY,
      version: DRAFT_VERSION,
      storage: createJSONStorage(() => localStorage),
      partialize: (state) => zDraft.parse(state),
      merge: (persisted, current) => {
        const stored = zStoredValues.safeParse(persisted);

        return stored.success ? { ...current, ...restoredDraft(stored.data, zDraft.parse(current)) } : current;
      },
    },
  ),
);
