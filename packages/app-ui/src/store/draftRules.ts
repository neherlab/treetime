import type { JobId, RunRecord, SparseConfig, UiDraft } from "@neherlab/app-contracts";

import { commandSettings } from "../settings/catalog";
import { normalizeConfig } from "../settings/config";
import { sameJson } from "../settings/json";
import { rerunDraft } from "../settings/rerun";
import { freshDraft } from "./draftSchema";

export function ownDraft(state: { draft: UiDraft; parked: FormSlot | undefined }): UiDraft {
  return state.parked?.draft ?? state.draft;
}

export function isPristine(draft: UiDraft): boolean {
  const fresh = freshDraft(draft.command);

  return (
    Object.keys(draft.sources).length === 0 &&
    draft.upload_run_id === undefined &&
    draft.from_run_id === undefined &&
    (draft.choices ?? []).length === 0 &&
    sameConfig(draft, fresh.config)
  );
}

export function sameConfig(draft: UiDraft, config: SparseConfig): boolean {
  const specs = commandSettings(draft.command).settings;

  return sameJson(normalizeConfig(specs, draft.config), normalizeConfig(specs, config));
}

export function runOptionsDraft(record: RunRecord, settings: SparseConfig, base: UiDraft): UiDraft | undefined {
  if (!record.prepared) {
    return undefined;
  }

  const { config, inputLabels } = rerunDraft(record, settings);

  const draft: UiDraft = {
    ...displayOf(base),
    command: record.command,
    config,
    sources: Object.fromEntries(
      Object.entries(inputLabels).map(([key, label]) => [key, { label, origin: { kind: "run" }, run_id: record.id }]),
    ),
    choices: [],
    from_run_id: record.id,
  };

  return record.example === undefined ? draft : { ...draft, example: record.example };
}

export function blankDraft(base: UiDraft, command: UiDraft["command"]): UiDraft {
  return { ...freshDraft(command), ...displayOf(base) };
}

export function withoutAbandonedUpload(draft: UiDraft, pendingPicks: number): UiDraft {
  const runId = draft.upload_run_id;

  if (
    runId === undefined ||
    pendingPicks > 0 ||
    Object.values(draft.sources).some((source) => source.run_id === runId)
  ) {
    return draft;
  }

  const { upload_run_id: _abandoned, ...kept } = draft;

  return kept;
}

export function lostOnDiscard(state: SlotState): boolean {
  const ownEdited = state.viewing === undefined && !isPristine(state.draft);

  return ownEdited || (state.parked !== undefined && !isPristine(state.parked.draft));
}

export function claimedSlot(state: SlotState, selection: symbol): "current" | "parked" {
  const parked = state.parked?.selections ?? {};

  return Object.values(parked).includes(selection) ? "parked" : "current";
}

export function slotHoldingKeys(state: SlotState, keys: readonly string[]): "current" | "parked" {
  const parked = state.parked?.selections ?? {};
  const [first] = keys;
  const selection = first === undefined ? undefined : parked[first];

  return selection !== undefined && keys.every((key) => parked[key] === selection) ? "parked" : "current";
}

export function heldUploadRuns(state: HeldDrafts): Set<JobId> {
  const drafts = [state.draft, state.parked?.draft, state.undo?.previous.draft, state.undo?.previous.parked].filter(
    (draft) => draft !== undefined,
  );

  return new Set(drafts.flatMap((draft) => (draft.upload_run_id === undefined ? [] : [draft.upload_run_id])));
}

export function releasedUploadRuns(before: HeldDrafts, after: HeldDrafts): JobId[] {
  const kept = heldUploadRuns(after);

  return [...heldUploadRuns(before)].filter((runId) => !kept.has(runId));
}

export function heldCreations(state: HeldCreations): Set<Promise<JobId>> {
  const creations = [state.uploadRunCreation, state.parked?.uploadRunCreation];

  return new Set(creations.filter((creation) => creation !== undefined));
}

export function restoredSnapshot(previous: DraftSnapshot, shownRunId: JobId | undefined): DraftSnapshot {
  if (previous.viewing === undefined || previous.viewing.runId === shownRunId) {
    return previous;
  }

  return { draft: previous.parked ?? previous.draft, viewing: undefined, parked: undefined };
}

export interface FormSlot {
  draft: UiDraft;
  selections: Readonly<Record<string, symbol>>;
  uploadRunCreation: Promise<JobId> | undefined;
}

export interface SlotState extends FormSlot {
  viewing: RunView | undefined;
  parked: FormSlot | undefined;
}

export interface DraftSnapshot {
  draft: UiDraft;
  viewing: RunView | undefined;
  parked: UiDraft | undefined;
}

export interface RunView {
  runId: JobId;
}

export interface HeldDrafts {
  draft: UiDraft;
  parked: FormSlot | undefined;
  undo: DraftUndo | undefined;
}

export interface HeldCreations {
  uploadRunCreation: Promise<JobId> | undefined;
  parked: FormSlot | undefined;
}

export interface DraftUndo {
  token: symbol;
  reason: DiscardReason;
  previous: DraftSnapshot;
}

export type DiscardReason = "reset" | "new" | "replaced" | "revert";

function displayOf({ view, search, changed_only, code_format }: UiDraft) {
  return { view, search, changed_only, code_format };
}
