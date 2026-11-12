import type { JobId, RunRecord, SparseConfig, UiDraft } from "@neherlab/app-contracts";

import { commandSettings } from "../settings/catalog";
import { normalizeConfig } from "../settings/config";
import { sameJson } from "../settings/json";
import { rerunDraft } from "../settings/rerun";
import { freshDraft } from "./draftSchema";

export function ownDraft(snapshot: DraftSnapshot): UiDraft {
  return snapshot.viewing?.parked ?? snapshot.draft;
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

export function runOptionsDraft(record: RunRecord, settings: SparseConfig, base: UiDraft): UiDraft {
  const { config, inputLabels } = rerunDraft(record, settings);

  return {
    ...displayOf(base),
    command: record.command,
    config,
    sources: Object.fromEntries(
      Object.entries(inputLabels).map(([key, label]) => [key, { label, origin: { kind: "run" }, run_id: record.id }]),
    ),
    choices: [],
    from_run_id: record.id,
  };
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

export function heldUploadRuns(state: HeldDrafts): Set<JobId> {
  const snapshots = state.undo === undefined ? [state] : [state, state.undo.previous];
  const drafts = snapshots.flatMap((snapshot) => [snapshot.draft, ownDraft(snapshot)]);

  return new Set(drafts.flatMap((draft) => (draft.upload_run_id === undefined ? [] : [draft.upload_run_id])));
}

export function releasedUploadRuns(before: HeldDrafts, after: HeldDrafts): JobId[] {
  const kept = heldUploadRuns(after);

  return [...heldUploadRuns(before)].filter((runId) => !kept.has(runId));
}

export function restoredSnapshot(previous: DraftSnapshot, shownRunId: JobId | undefined): DraftSnapshot {
  if (previous.viewing === undefined || previous.viewing.runId === shownRunId) {
    return previous;
  }

  return { draft: previous.viewing.parked, viewing: undefined };
}

export interface DraftSnapshot {
  draft: UiDraft;
  viewing: RunView | undefined;
}

export interface RunView {
  runId: JobId;
  parked: UiDraft;
}

export interface HeldDrafts extends DraftSnapshot {
  undo: DraftUndo | undefined;
}

export interface DraftUndo {
  token: symbol;
  reason: DiscardReason;
  previous: DraftSnapshot;
}

export type DiscardReason = "reset" | "new" | "replaced";

function displayOf({ view, search, changed_only, code_format }: UiDraft) {
  return { view, search, changed_only, code_format };
}
