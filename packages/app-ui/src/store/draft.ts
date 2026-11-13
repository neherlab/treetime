import type { ActiveChoice, AppCommand, JobId, SparseConfig, UiDraft, UiDraftSource } from "@neherlab/app-contracts";
import { create } from "zustand";

import { withDraftExample } from "./draftExample";
import {
  blankDraft,
  type DiscardReason,
  type DraftSnapshot,
  type DraftUndo,
  type FormSlot,
  isPristine,
  lostOnDiscard,
  restoredSnapshot,
  type RunView,
  sameConfig,
  slotHoldingKeys,
  withoutAbandonedUpload,
} from "./draftRules";
import { forgetRuns, freshDraft, type LostInputs, type QueuedExample, withPickedChoice } from "./draftSchema";

interface DraftState extends FormSlot {
  viewing: RunView | undefined;
  parked: FormSlot | undefined;
  lastShown: JobId | undefined;
  undo: DraftUndo | undefined;
  queuedExamples: readonly QueuedExample[];
  epoch: number;
  busyRuns: readonly JobId[];
  focus: SettingFocus | undefined;
}

export interface SettingFocus {
  key: string;
  token: symbol;
}

interface DraftActions {
  load: (draft: UiDraft) => void;
  update: (changes: Partial<Pick<UiDraft, "view" | "search" | "changed_only" | "code_format">>) => void;
  setConfig: (config: SparseConfig) => void;
  pickChoice: (picked: ActiveChoice) => void;
  setSource: (key: string, source: UiDraftSource | undefined) => void;
  assignInputs: (changes: Pick<UiDraft, "config" | "sources">, keys: readonly string[]) => void;
  claimSelections: (keys: readonly string[], selection: symbol) => void;
  beginUploadRun: (creation: Promise<JobId>, selection: symbol | undefined) => void;
  settleUploadRun: (creation: Promise<JobId>, id: JobId | undefined) => void;
  beginRunWork: (runId: JobId) => void;
  endRunWork: (runId: JobId) => void;
  finishStart: () => void;
  viewRun: (runId: JobId, options: UiDraft) => void;
  refreshView: (runId: JobId, options: UiDraft) => void;
  leaveRun: () => void;
  discard: (reason: Exclude<DiscardReason, "replaced" | "revert">, command: AppCommand) => void;
  revertToRun: (runId: JobId, options: UiDraft) => void;
  forgetInputs: (lost: LostInputs) => void;
  focusSetting: (key: string) => void;
  clearFocus: (token: symbol) => void;
  undoDiscard: (token: symbol, shownRunId: JobId | undefined) => void;
  dropUndo: (token: symbol) => void;
  queueExamples: (examples: readonly QueuedExample[]) => void;
  takeQueuedExamples: () => void;
}

type Store = DraftState & DraftActions;

export const useDraftStore = create<Store>()((set) => ({
  ...emptySlot(freshDraft("timetree")),
  viewing: undefined,
  parked: undefined,
  lastShown: undefined,
  undo: undefined,
  queuedExamples: [],
  epoch: 0,
  busyRuns: [],
  focus: undefined,
  load: (draft) => {
    set((state) => ({
      ...edited(withoutAbandonedUpload(withDraftExample(draft, draft.example), 0)),
      epoch: state.epoch + 1,
      selections: {},
      uploadRunCreation: undefined,
    }));
  },
  update: (changes) => {
    set(({ draft }) => ({ draft: { ...draft, ...changes } }));
  },
  setConfig: (config) => {
    set((state) => (sameConfig(state.draft, config) ? {} : edited({ ...state.draft, config })));
  },
  pickChoice: (picked) => {
    set((state) => edited({ ...state.draft, choices: withPickedChoice(state.draft.choices, picked) }));
  },
  setSource: (key, source) => {
    set((state) => {
      const sources = Object.fromEntries(Object.entries(state.draft.sources).filter(([name]) => name !== key));
      const selections = Object.fromEntries(Object.entries(state.selections).filter(([name]) => name !== key));
      const next = { ...state.draft, sources: source === undefined ? sources : { ...sources, [key]: source } };
      const draft = withoutAbandonedUpload(withDraftExample(next, state.draft.example), Object.keys(selections).length);

      return { ...edited(draft), selections };
    });
  },
  assignInputs: (changes, keys) => {
    set((state) => {
      const target = slotHoldingKeys(state, keys) === "parked" ? state.parked : undefined;
      const slot = target ?? state;
      const selections = Object.fromEntries(Object.entries(slot.selections).filter(([name]) => !keys.includes(name)));
      const next = withDraftExample({ ...slot.draft, ...changes }, slot.draft.example);
      const draft = withoutAbandonedUpload(next, Object.keys(selections).length);

      if (target !== undefined) {
        return { parked: { ...target, draft, selections } };
      }

      return { ...edited(draft), epoch: state.epoch + 1, selections };
    });
  },
  claimSelections: (keys, selection) => {
    set(({ selections }) => ({
      selections: { ...selections, ...Object.fromEntries(keys.map((key) => [key, selection])) },
    }));
  },
  beginUploadRun: (creation, selection) => {
    set((state) => {
      const parked = state.parked;

      if (parked !== undefined && selection !== undefined && Object.values(parked.selections).includes(selection)) {
        return { parked: { ...parked, uploadRunCreation: creation } };
      }

      return { uploadRunCreation: creation };
    });
  },
  settleUploadRun: (creation, id) => {
    set((state) => {
      const settle = (slot: FormSlot): FormSlot => ({
        ...slot,
        uploadRunCreation: undefined,
        draft: id === undefined ? slot.draft : { ...slot.draft, upload_run_id: id },
      });

      if (state.uploadRunCreation === creation) {
        return settle(state);
      }

      return state.parked?.uploadRunCreation === creation ? { parked: settle(state.parked) } : {};
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
      ...(state.parked ?? emptySlot(blankDraft(state.draft, state.draft.command))),
      viewing: undefined,
      parked: undefined,
      epoch: state.epoch + 1,
    }));
  },
  viewRun: (runId, options) => {
    set((state) => {
      if (state.viewing?.runId === runId || (state.viewing === undefined && state.draft.from_run_id === runId)) {
        return { lastShown: runId };
      }

      const own = state.viewing === undefined ? slotOf(state) : undefined;
      const replaced = state.parked !== undefined && own !== undefined && !isPristine(own.draft);

      return {
        ...emptySlot(options),
        viewing: { runId },
        parked: state.parked ?? own,
        lastShown: runId,
        undo: replaced ? undoOf("replaced", state) : state.undo,
        epoch: state.epoch + 1,
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
            ...(state.parked ?? emptySlot(blankDraft(state.draft, state.draft.command))),
            viewing: undefined,
            parked: undefined,
            lastShown: undefined,
            epoch: state.epoch + 1,
          },
    );
  },
  discard: (reason, command) => {
    set((state) => {
      if (reason === "reset" && state.parked !== undefined) {
        const variantLost = state.viewing === undefined && !isPristine(state.draft);

        return {
          ...state.parked,
          viewing: undefined,
          parked: undefined,
          undo: variantLost ? undoOf(reason, state) : state.undo,
          epoch: state.epoch + 1,
        };
      }

      return {
        ...emptySlot(blankDraft(state.draft, command)),
        viewing: undefined,
        parked: undefined,
        undo: lostOnDiscard(state) ? undoOf(reason, state) : state.undo,
        epoch: state.epoch + 1,
      };
    });
  },
  revertToRun: (runId, options) => {
    set((state) => ({
      ...emptySlot(options),
      viewing: { runId },
      parked: state.parked ?? emptySlot(blankDraft(state.draft, state.draft.command)),
      lastShown: runId,
      undo: undoOf("revert", state),
      epoch: state.epoch + 1,
    }));
  },
  forgetInputs: (lost) => {
    set((state) => {
      const current = forgetRuns(state.draft, lost).draft;

      const changes: Partial<DraftState> = sameDraft(current, state.draft)
        ? {}
        : { draft: current, epoch: state.epoch + 1 };

      if (state.parked === undefined) {
        return changes;
      }

      const parked = forgetRuns(state.parked.draft, lost).draft;

      return sameDraft(parked, state.parked.draft)
        ? changes
        : { ...changes, parked: { ...state.parked, draft: parked } };
    });
  },
  focusSetting: (key) => {
    set({ focus: { key, token: Symbol(key) } });
  },
  clearFocus: (token) => {
    set((state) => (state.focus?.token === token ? { focus: undefined } : {}));
  },
  undoDiscard: (token, shownRunId) => {
    set((state) => {
      if (state.undo?.token !== token) {
        return {};
      }

      const restored = restoredSnapshot(state.undo.previous, shownRunId);

      return {
        ...emptySlot(restored.draft),
        viewing: restored.viewing,
        parked: restored.parked === undefined ? undefined : emptySlot(restored.parked),
        undo: undefined,
        epoch: state.epoch + 1,
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

function edited(draft: UiDraft): Pick<DraftState, "draft" | "viewing"> {
  return { draft, viewing: undefined };
}

function emptySlot(draft: UiDraft): FormSlot {
  return { draft, selections: {}, uploadRunCreation: undefined };
}

function slotOf({ draft, selections, uploadRunCreation }: FormSlot): FormSlot {
  return { draft, selections, uploadRunCreation };
}

function undoOf(reason: DiscardReason, state: DraftState): DraftUndo {
  return { token: Symbol(reason), reason, previous: snapshotOf(state) };
}

function snapshotOf({ draft, viewing, parked }: DraftState): DraftSnapshot {
  return { draft, viewing, parked: parked?.draft };
}

function sameDraft(left: UiDraft, right: UiDraft): boolean {
  return JSON.stringify(left) === JSON.stringify(right);
}
