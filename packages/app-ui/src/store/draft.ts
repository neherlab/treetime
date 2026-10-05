import type { AppCommand, JobId, SparseConfig, UiDraft, UiDraftSource } from "@neherlab/app-contracts";
import { create } from "zustand";

import { freshDraft } from "./draftSchema";

interface DraftState {
  draft: UiDraft;
  epoch: number;
}

interface DraftActions {
  load: (draft: UiDraft) => void;
  update: (changes: Partial<Pick<UiDraft, "view" | "search" | "changed_only" | "code_format">>) => void;
  setConfig: (config: SparseConfig) => void;
  setSource: (key: string, source: UiDraftSource | undefined) => void;
  setUploadRun: (id: JobId) => void;
  clearRuns: () => void;
  reset: (command: AppCommand) => void;
}

export const useDraftStore = create<DraftState & DraftActions>()((set) => ({
  draft: freshDraft("timetree"),
  epoch: 0,
  load: (draft) => {
    set((state) => ({ draft, epoch: state.epoch + 1 }));
  },
  update: (changes) => {
    set(({ draft }) => ({ draft: { ...draft, ...changes } }));
  },
  setConfig: (config) => {
    set(({ draft }) => ({ draft: { ...draft, config } }));
  },
  setSource: (key, source) => {
    set(({ draft }) => {
      const sources = Object.fromEntries(Object.entries(draft.sources).filter(([name]) => name !== key));

      return { draft: { ...draft, sources: source === undefined ? sources : { ...sources, [key]: source } } };
    });
  },
  setUploadRun: (id) => {
    set(({ draft }) => ({ draft: { ...draft, upload_run_id: id } }));
  },
  clearRuns: () => {
    set(({ draft: { from_run_id: _fromRun, upload_run_id: _uploadRun, ...draft } }) => ({ draft }));
  },
  reset: (command) => {
    set((state) => ({ draft: freshDraft(command), epoch: state.epoch + 1 }));
  },
}));
