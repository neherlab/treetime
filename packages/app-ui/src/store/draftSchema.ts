import { zUiDraft, type AppCommand, type JobId, type UiDraft } from "@neherlab/app-contracts";
import type * as z from "zod";

import { COMMAND_SETTINGS } from "../settings/catalog";
import { defaultConfig } from "../settings/config";
import { zJsonObject, type JsonObject } from "../settings/json";

export type Draft = Omit<UiDraft, "config" | "from_run_id" | "upload_run_id"> & {
  config: JsonObject;
  from_run_id: JobId | null;
  upload_run_id: JobId | null;
};

export function freshDraft(command: AppCommand): Draft {
  return {
    command,
    config: defaultConfig(COMMAND_SETTINGS[command].specs),
    sources: {},
    from_run_id: null,
    upload_run_id: null,
    view: "main",
    search: "",
    changed_only: false,
    code_format: "cli",
  };
}

export function storedDraft(draft: z.output<typeof zUiDraft>): Draft {
  return {
    ...draft,
    config: zJsonObject.parse(draft.config),
    sources: Object.fromEntries(
      Object.entries(draft.sources).map(([key, { label, origin, size }]) => [
        key,
        { label, origin, size: size ?? null },
      ]),
    ),
    from_run_id: draft.from_run_id ?? null,
    upload_run_id: draft.upload_run_id ?? null,
  };
}
