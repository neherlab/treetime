import type { ActiveChoice, AppCommand, JobId, UiDraft } from "@neherlab/app-contracts";

import { commandSettings } from "../settings/catalog";
import { defaultConfig, resetValue } from "../settings/config";
import { isListInput, pathList } from "../settings/inputs";
import { withSetting } from "../settings/json";

export function freshDraft(command: AppCommand): UiDraft {
  return {
    command,
    config: defaultConfig(commandSettings(command).settings),
    sources: {},
    choices: [],
    view: "main",
    search: "",
    changed_only: false,
    code_format: "cli",
  };
}

export function withPickedChoice(choices: readonly ActiveChoice[] | undefined, picked: ActiveChoice): ActiveChoice[] {
  return [...(choices ?? []).filter((active) => active.choice !== picked.choice), picked];
}

export function referencedRuns(draft: UiDraft): JobId[] {
  const ids = [draft.upload_run_id, draft.from_run_id, ...Object.values(draft.sources).map((source) => source.run_id)];

  return [...new Set(ids.filter((id) => id !== undefined))];
}

export function forgetRuns(draft: UiDraft, lost: LostInputs): ForgottenRuns {
  const dropped = Object.entries(draft.sources).filter(([key]) =>
    pathList(draft.config[key]).some((path) => lost.paths.has(path)),
  );

  const keys = new Set(dropped.map(([key]) => key));

  const config = commandSettings(draft.command)
    .settings.filter((spec) => keys.has(spec.key))
    .reduce((current, spec) => withSetting(current, spec.key, resetValue(spec)), draft.config);

  const sources = Object.fromEntries(Object.entries(draft.sources).filter(([key]) => !keys.has(key)));
  const { upload_run_id: uploadRun, from_run_id: fromRun, ...kept } = draft;
  const forgotten: UiDraft = { ...kept, config, sources };

  if (uploadRun !== undefined && !lost.runs.has(uploadRun)) {
    forgotten.upload_run_id = uploadRun;
  }

  if (fromRun !== undefined && !lost.runs.has(fromRun)) {
    forgotten.from_run_id = fromRun;
  }

  return {
    draft: forgotten,
    dropped: dropped.flatMap(([, source]) => (source.origin.kind === "example" ? [] : [source.label])),
    examples: dropped.flatMap(([key, { origin }]) =>
      origin.kind === "example" ? [{ key, list: isListInput(draft.command, key), path: origin.path }] : [],
    ),
  };
}

export function draftInputPaths(draft: UiDraft): string[] {
  return [...new Set(Object.keys(draft.sources).flatMap((key) => pathList(draft.config[key])))];
}

export function runSourcePaths(draft: UiDraft, runs: ReadonlySet<JobId>): string[] {
  return Object.entries(draft.sources).flatMap(([key, source]) =>
    source.run_id !== undefined && runs.has(source.run_id) ? pathList(draft.config[key]) : [],
  );
}

export interface LostInputs {
  runs: ReadonlySet<JobId>;
  paths: ReadonlySet<string>;
}

export interface ForgottenRuns {
  draft: UiDraft;
  dropped: string[];
  examples: QueuedExample[];
}

export interface QueuedExample {
  key: string;
  list: boolean;
  path: string;
}
