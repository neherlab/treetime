import type { AppCommand } from "@neherlab/app-contracts";
import * as z from "zod";
import { create } from "zustand";
import { createJSONStorage, persist } from "zustand/middleware";

import { COMMAND_SETTINGS } from "../settings/catalog";
import { APP_COMMANDS } from "../settings/commands";
import { defaultConfig } from "../settings/config";
import { zJsonObject, type JsonObject } from "../settings/json";

export type InputOrigin = "dataset" | "upload" | "local" | "run" | "config";

export interface InputSource {
  label: string;
  origin: InputOrigin;
  size: number | null;
}

export type SettingsView = "main" | "all";

export type CodeFormat = "cli" | "yaml";

interface DraftState {
  command: AppCommand;
  config: JsonObject;
  sources: Record<string, InputSource>;
  title: string;
  fromRunId: string | null;
  uploadRunId: string | null;
  view: SettingsView;
  search: string;
  changedOnly: boolean;
  codeFormat: CodeFormat;
  epoch: number;
}

interface DraftActions {
  load: (draft: Partial<Omit<DraftState, "epoch">>) => void;
  update: (draft: Partial<Omit<DraftState, "epoch" | "command" | "config">>) => void;
  setConfig: (config: JsonObject) => void;
  setSource: (key: string, source: InputSource | null) => void;
  reset: (command: AppCommand) => void;
}

const DRAFT_STORAGE_KEY = "treetime-draft";

const DRAFT_VERSION = 1;

const zInputSource = z.object({
  label: z.string(),
  origin: z.enum(["dataset", "upload", "local", "run", "config"]),
  size: z.number().nullable(),
});

const zPersistedDraft = z
  .object({
    command: z.enum(APP_COMMANDS),
    config: zJsonObject,
    sources: z.record(z.string(), zInputSource),
    title: z.string(),
    fromRunId: z.string().nullable(),
    uploadRunId: z.string().nullable(),
    view: z.enum(["main", "all"]),
    search: z.string(),
    changedOnly: z.boolean(),
    codeFormat: z.enum(["cli", "yaml"]),
  })
  .partial();

function freshDraft(command: AppCommand): Omit<DraftState, "epoch"> {
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

export const useDraftStore = create<DraftState & DraftActions>()(
  persist(
    (set) => ({
      ...freshDraft("timetree"),
      epoch: 0,
      load: (draft) => {
        set((state) => ({ ...draft, epoch: state.epoch + 1 }));
      },
      update: (draft) => {
        set(draft);
      },
      setConfig: (config) => {
        set({ config });
      },
      setSource: (key, source) => {
        set((state) => {
          const sources = Object.fromEntries(Object.entries(state.sources).filter(([name]) => name !== key));

          return { sources: source === null ? sources : { ...sources, [key]: source } };
        });
      },
      reset: (command) => {
        set((state) => ({ ...freshDraft(command), epoch: state.epoch + 1 }));
      },
    }),
    {
      name: DRAFT_STORAGE_KEY,
      version: DRAFT_VERSION,
      storage: createJSONStorage(() => localStorage),
      partialize: (state) => ({
        command: state.command,
        config: state.config,
        sources: state.sources,
        title: state.title,
        fromRunId: state.fromRunId,
        uploadRunId: state.uploadRunId,
        view: state.view,
        search: state.search,
        changedOnly: state.changedOnly,
        codeFormat: state.codeFormat,
      }),
      merge: (persisted, current) => {
        const parsed = zPersistedDraft.safeParse(persisted);

        if (!parsed.success) {
          return current;
        }

        const draft = parsed.data;

        return {
          ...current,
          command: draft.command ?? current.command,
          config: draft.config ?? current.config,
          sources: draft.sources ?? current.sources,
          title: draft.title ?? current.title,
          fromRunId: draft.fromRunId ?? current.fromRunId,
          uploadRunId: draft.uploadRunId ?? current.uploadRunId,
          view: draft.view ?? current.view,
          search: draft.search ?? current.search,
          changedOnly: draft.changedOnly ?? current.changedOnly,
          codeFormat: draft.codeFormat ?? current.codeFormat,
        };
      },
    },
  ),
);
