import type { AppCommand, JobId, SparseConfig, UiDraft, UiDraftSource } from "@neherlab/app-contracts";
import { create } from "zustand";

import { freshDraft } from "./draftSchema";

interface DraftState {
  draft: UiDraft;
  epoch: number;
  selections: Readonly<Record<string, symbol>>;
  uploadRunCreation: Promise<JobId> | undefined;
}

interface DraftActions {
  load: (draft: UiDraft) => void;
  update: (changes: Partial<Pick<UiDraft, "view" | "search" | "changed_only" | "code_format">>) => void;
  setConfig: (config: SparseConfig) => void;
  setSource: (key: string, source: UiDraftSource | undefined) => void;
  claimSelections: (keys: readonly string[], selection: symbol) => void;
  beginUploadRun: (creation: Promise<JobId>) => void;
  settleUploadRun: (creation: Promise<JobId>, id: JobId | undefined) => void;
  clearRuns: () => void;
  reset: (command: AppCommand) => void;
}

export const useDraftStore = create<DraftState & DraftActions>()((set) => ({
  draft: freshDraft("timetree"),
  epoch: 0,
  selections: {},
  uploadRunCreation: undefined,
  load: (draft) => {
    set((state) => ({ draft, epoch: state.epoch + 1, selections: {}, uploadRunCreation: undefined }));
  },
  update: (changes) => {
    set(({ draft }) => ({ draft: { ...draft, ...changes } }));
  },
  setConfig: (config) => {
    set(({ draft }) => ({ draft: { ...draft, config } }));
  },
  setSource: (key, source) => {
    set(({ draft, selections }) => {
      const sources = Object.fromEntries(Object.entries(draft.sources).filter(([name]) => name !== key));

      return {
        draft: { ...draft, sources: source === undefined ? sources : { ...sources, [key]: source } },
        selections: Object.fromEntries(Object.entries(selections).filter(([name]) => name !== key)),
      };
    });
  },
  claimSelections: (keys, selection) => {
    set(({ selections }) => ({
      selections: { ...selections, ...Object.fromEntries(keys.map((key) => [key, selection])) },
    }));
  },
  beginUploadRun: (creation) => {
    set({ uploadRunCreation: creation });
  },
  settleUploadRun: (creation, id) => {
    set((state) => {
      if (state.uploadRunCreation !== creation) {
        return {};
      }

      return {
        uploadRunCreation: undefined,
        draft: id === undefined ? state.draft : { ...state.draft, upload_run_id: id },
      };
    });
  },
  clearRuns: () => {
    set(({ draft: { from_run_id: _fromRun, upload_run_id: _uploadRun, ...draft } }) => ({
      draft,
      selections: {},
      uploadRunCreation: undefined,
    }));
  },
  reset: (command) => {
    set((state) => ({
      draft: freshDraft(command),
      epoch: state.epoch + 1,
      selections: {},
      uploadRunCreation: undefined,
    }));
  },
}));
