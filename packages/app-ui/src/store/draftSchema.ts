import { zAppCommand, type AppCommand } from "@neherlab/app-contracts";
import * as z from "zod";

import { COMMAND_SETTINGS } from "../settings/catalog";
import { defaultConfig } from "../settings/config";
import { zJsonObject } from "../settings/json";

const zInputOrigin = z.enum(["dataset", "upload", "local", "run", "config"]);

const zInputSource = z.object({
  label: z.string(),
  origin: zInputOrigin,
  size: z.number().nullable(),
});

const zSettingsView = z.enum(["main", "all"]);

const zCodeFormat = z.enum(["cli", "yaml"]);

export const zDraft = z.object({
  command: zAppCommand,
  config: zJsonObject,
  sources: z.record(z.string(), zInputSource),
  title: z.string(),
  fromRunId: z.string().nullable(),
  uploadRunId: z.string().nullable(),
  view: zSettingsView,
  search: z.string(),
  changedOnly: z.boolean(),
  codeFormat: zCodeFormat,
});

export const zStoredValues = z.record(z.string(), z.unknown());

export type InputOrigin = z.infer<typeof zInputOrigin>;

export type InputSource = z.infer<typeof zInputSource>;

export type SettingsView = z.infer<typeof zSettingsView>;

export type CodeFormat = z.infer<typeof zCodeFormat>;

export type Draft = z.infer<typeof zDraft>;

type StoredValues = z.infer<typeof zStoredValues>;

export function restoredDraft(stored: StoredValues, current: Draft): Draft {
  const parsed = zDraft.safeParse({ ...current, ...stored });

  return parsed.success ? parsed.data : current;
}

export function freshDraft(command: AppCommand): Draft {
  return {
    command,
    config: defaultConfig(COMMAND_SETTINGS[command].specs),
    sources: {},
    title: "",
    fromRunId: null,
    uploadRunId: null,
    view: "main",
    search: "",
    changedOnly: false,
    codeFormat: "cli",
  };
}
