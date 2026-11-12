import type { ActiveChoice, AppCommand, JobId, SparseConfig, UiDraft, UiDraftSource } from "@neherlab/app-contracts";
import { create } from "zustand";

import {
  blankDraft,
  type DiscardReason,
  type DraftUndo,
  isPristine,
  ownDraft,
  restoredSnapshot,
  type RunView,
  sameConfig,
  withoutAbandonedUpload,
} from "./draftRules";
import { freshDraft, type QueuedExample, withPickedChoice } from "./draftSchema";

interface DraftState {
  draft: UiDraft;
  viewing: RunView | undefined;
  lastShown: JobId | undefined;
  undo: DraftUndo | undefined;
  queuedExamples: readonly QueuedExample[];
  epoch: number;
  selections: Readonly<Record<string, symbol>>;
  uploadRunCreation: Promise<JobId> | undefined;
  busyRuns: readonly JobId[];
}

interface DraftActions {
  load: (draft: UiDraft) => void;
  update: (changes: Partial<Pick<UiDraft, "view" | "search" | "changed_only" | "code_format">>) => void;
  setConfig: (config: SparseConfig) => void;
  pickChoice: (picked: ActiveChoice) => void;
  setSource: (key: string, source: UiDraftSource | undefined) => void;
  assignInputs: (changes: Pick<UiDraft, "config" | "sources">, keys: readonly string[]) => void;
  claimSelections: (keys: readonly string[], selection: symbol) => void;
  beginUploadRun: (creation: Promise<JobId>) => void;
  settleUploadRun: (creation: Promise<JobId>, id: JobId | undefined) => void;
  beginRunWork: (runId: JobId) => void;
  endRunWork: (runId: JobId) => void;
  finishStart: () => void;
  viewRun: (runId: JobId, options: UiDraft) => void;
  refreshView: (runId: JobId, options: UiDraft) => void;
  leaveRun: () => void;
  discard: (reason: Exclude<DiscardReason, "replaced">, command: AppCommand) => void;
  undoDiscard: (token: symbol, shownRunId: JobId | undefined) => void;
  dropUndo: (token: symbol) => void;
  queueExamples: (examples: readonly QueuedExample[]) => void;
  takeQueuedExamples: () => void;
}

type Store = DraftState & DraftActions;

export const useDraftStore = create<Store>()((set) => ({
  draft: freshDraft("timetree"),
  viewing: undefined,
  lastShown: undefined,
  undo: undefined,
  queuedExamples: [],
  epoch: 0,
  selections: {},
  uploadRunCreation: undefined,
  busyRuns: [],
  load: (draft) => {
    set((state) => ({
      ...edited(state, withoutAbandonedUpload(draft, 0)),
      epoch: state.epoch + 1,
      selections: {},
      uploadRunCreation: undefined,
    }));
  },
  update: (changes) => {
    set(({ draft }) => ({ draft: { ...draft, ...changes } }));
  },
  setConfig: (config) => {
    set((state) => (sameConfig(state.draft, config) ? {} : edited(state, { ...state.draft, config })));
  },
  pickChoice: (picked) => {
    set((state) => edited(state, { ...state.draft, choices: withPickedChoice(state.draft.choices, picked) }));
  },
  setSource: (key, source) => {
    set((state) => {
      const sources = Object.fromEntries(Object.entries(state.draft.sources).filter(([name]) => name !== key));
      const selections = Object.fromEntries(Object.entries(state.selections).filter(([name]) => name !== key));

      const draft = withoutAbandonedUpload(
        { ...state.draft, sources: source === undefined ? sources : { ...sources, [key]: source } },
        Object.keys(selections).length,
      );

      return { ...edited(state, draft), selections };
    });
  },
  assignInputs: (changes, keys) => {
    set((state) => {
      const selections = Object.fromEntries(Object.entries(state.selections).filter(([name]) => !keys.includes(name)));
      const draft = withoutAbandonedUpload({ ...state.draft, ...changes }, Object.keys(selections).length);

      return { ...edited(state, draft), epoch: state.epoch + 1, selections };
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
  beginRunWork: (runId) => {
    set(({ busyRuns }) => ({ busyRuns: [...busyRuns, runId] }));
  },
  endRunWork: (runId) => {
    set(({ busyRuns }) => {
      const index = busyRuns.indexOf(runId);

      return index < 0 ? {} : { busyRuns: busyRuns.toSpliced(index, 1) };
    });
  },
  finishStart: () => {
    set((state) => ({
      draft: blankDraft(state.draft, state.draft.command),
      viewing: undefined,
      epoch: state.epoch + 1,
      selections: {},
      uploadRunCreation: undefined,
    }));
  },
  viewRun: (runId, options) => {
    set((state) => {
      if (state.viewing?.runId === runId || (state.viewing === undefined && state.draft.from_run_id === runId)) {
        return { lastShown: runId };
      }

      return {
        draft: options,
        viewing: { runId, parked: ownDraft(state) },
        lastShown: runId,
        epoch: state.epoch + 1,
        selections: {},
        uploadRunCreation: undefined,
      };
    });
  },
  refreshView: (runId, options) => {
    set((state) => {
      if (state.viewing?.runId !== runId || sameDraft(state.draft, options)) {
        return {};
      }

      return { draft: options, epoch: state.epoch + 1 };
    });
  },
  leaveRun: () => {
    set((state) =>
      state.viewing === undefined
        ? { lastShown: undefined }
        : {
            draft: state.viewing.parked,
            viewing: undefined,
            lastShown: undefined,
            epoch: state.epoch + 1,
            selections: {},
            uploadRunCreation: undefined,
          },
    );
  },
  discard: (reason, command) => {
    set((state) => {
      const kept = isPristine(state.draft) && isPristine(ownDraft(state));

      return {
        draft: blankDraft(state.draft, command),
        viewing: undefined,
        undo: kept ? state.undo : { token: Symbol(reason), reason, previous: snapshotOf(state) },
        epoch: state.epoch + 1,
        selections: {},
        uploadRunCreation: undefined,
      };
    });
  },
  undoDiscard: (token, shownRunId) => {
    set((state) => {
      if (state.undo?.token !== token) {
        return {};
      }

      const restored = restoredSnapshot(state.undo.previous, shownRunId);

      return {
        ...restored,
        undo: undefined,
        epoch: state.epoch + 1,
        selections: {},
        uploadRunCreation: undefined,
      };
    });
  },
  dropUndo: (token) => {
    set((state) => (state.undo?.token === token ? { undo: undefined } : {}));
  },
  queueExamples: (examples) => {
    set({ queuedExamples: examples });
  },
  takeQueuedExamples: () => {
    set({ queuedExamples: [] });
  },
}));

function edited(state: DraftState, draft: UiDraft): Pick<DraftState, "draft" | "viewing" | "undo"> {
  if (state.viewing === undefined) {
    return { draft, viewing: undefined, undo: state.undo };
  }

  const replaced = !isPristine(state.viewing.parked);

  return {
    draft,
    viewing: undefined,
    undo: replaced ? { token: Symbol("replaced"), reason: "replaced", previous: snapshotOf(state) } : state.undo,
  };
}

function snapshotOf({ draft, viewing }: DraftState) {
  return { draft, viewing };
}

function sameDraft(left: UiDraft, right: UiDraft): boolean {
  return JSON.stringify(left) === JSON.stringify(right);
}
