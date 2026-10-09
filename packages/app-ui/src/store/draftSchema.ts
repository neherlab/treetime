import type { AppCommand, JobId, UiDraft } from "@neherlab/app-contracts";

import { commandSettings } from "../settings/catalog";
import { defaultConfig, resetValue } from "../settings/config";
import { setAt } from "../settings/json";

export function freshDraft(command: AppCommand): UiDraft {
  return {
    command,
    config: defaultConfig(commandSettings(command).settings),
    sources: {},
    view: "main",
    search: "",
    changed_only: false,
    code_format: "cli",
  };
}

export function referencedRuns(draft: UiDraft): JobId[] {
  const ids = [draft.upload_run_id, draft.from_run_id, ...Object.values(draft.sources).map((source) => source.run_id)];

  return [...new Set(ids.filter((id) => id !== undefined))];
}

export function forgetRuns(draft: UiDraft, lost: ReadonlySet<JobId>): ForgottenRuns {
  const dropped = Object.entries(draft.sources).filter(
    ([, source]) => source.run_id !== undefined && lost.has(source.run_id),
  );
  const keys = new Set(dropped.map(([key]) => key));
  const config = commandSettings(draft.command)
    .settings.filter((spec) => keys.has(spec.key))
    .reduce((current, spec) => setAt(current, spec.path, resetValue(spec)), draft.config);
  const sources = Object.fromEntries(Object.entries(draft.sources).filter(([key]) => !keys.has(key)));
  const { upload_run_id: uploadRun, from_run_id: fromRun, ...kept } = draft;

  return {
    draft: {
      ...kept,
      config,
      sources,
      ...(uploadRun === undefined || lost.has(uploadRun) ? {} : { upload_run_id: uploadRun }),
      ...(fromRun === undefined || lost.has(fromRun) ? {} : { from_run_id: fromRun }),
    },
    dropped: dropped.map(([, source]) => source.label),
  };
}

export interface ForgottenRuns {
  draft: UiDraft;
  dropped: string[];
}
