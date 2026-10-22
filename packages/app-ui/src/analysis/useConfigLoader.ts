import type { AppCommand, CheckConfigRequest, SparseConfig, UiDraftSource } from "@neherlab/app-contracts";
import { configCheck } from "@neherlab/app-contracts/client";
import { useNavigate } from "@tanstack/react-router";
import { useCallback } from "react";

import { useApiContext } from "../api/context";
import { commandSettings } from "../settings/catalog";
import { normalizeConfig, settingValue } from "../settings/config";
import { baseName, pathList } from "../settings/inputs";
import { sameJson } from "../settings/json";
import { useDraftStore } from "../store/draft";
import type { UploadedInputs } from "./exampleData";
import { releasedPaths } from "./fileSelection";
import { useRunUploads } from "./useRunUploads";

export type ConfigLoadResult = { loaded: true; command: AppCommand } | { loaded: false; messages: string[] };

export const SUPERSEDED: ConfigLoadResult = {
  loaded: false,
  messages: ["A later choice of inputs replaced this config"],
};

export interface ConfigLoadOptions {
  keepInputs?: boolean;
  folder?: string | undefined;
  uploaded?: UploadedInputs;
  isCurrent?: () => boolean;
}

export function useConfigLoader() {
  const { client } = useApiContext();
  const navigate = useNavigate();
  const { release } = useRunUploads();

  return useCallback(
    async (
      text: string,
      fallbackCommand: AppCommand,
      { keepInputs = false, folder, uploaded, isCurrent }: ConfigLoadOptions = {},
    ): Promise<ConfigLoadResult> => {
      const { draft, load } = useDraftStore.getState();
      const inputs = keepInputs ? inputSettings(draft.command, draft.config) : {};

      const body: CheckConfigRequest = { command: fallbackCommand, text, inputs };

      if (folder !== undefined) {
        body.folder = folder;
      }

      if (uploaded !== undefined) {
        body.replace_inputs = uploaded.replace_inputs;
      }

      const { data: result } = await configCheck({ client, body, throwOnError: true });

      if (result.status === "invalid") {
        return { loaded: false, messages: result.messages };
      }

      if (isCurrent?.() === false) {
        return SUPERSEDED;
      }

      const command = result.command;
      const specs = commandSettings(command).settings;
      const config = normalizeConfig(specs, result.config);
      const sources: Record<string, UiDraftSource> = {};
      const released: string[] = [];

      for (const spec of specs.filter((candidate) => candidate.role === "input")) {
        const value = settingValue(config, spec);
        const paths = pathList(value);
        const previous = draft.sources[spec.key];
        const upload = uploaded?.sources[spec.key];

        if (paths.length > 0) {
          sources[spec.key] =
            upload !== undefined && sameJson(value, uploaded?.replace_inputs[spec.key])
              ? upload
              : previousSource(previous, sameJson(value, draft.config[spec.key]), paths);
        }

        if (previous?.origin === "upload" && previous.run_id === draft.upload_run_id) {
          released.push(...releasedPaths(draft.config[spec.key], value));
        }
      }

      const { from_run_id: _fromRun, ...current } = draft;

      load({ ...current, command, config, sources });

      if (draft.upload_run_id !== undefined) {
        release(draft.upload_run_id, released);
      }

      await navigate({ to: "/new" });

      return { loaded: true, command };
    },
    [client, navigate, release],
  );
}

function previousSource(previous: UiDraftSource | undefined, unchanged: boolean, paths: string[]): UiDraftSource {
  return previous !== undefined && unchanged ? previous : { label: paths.map(baseName).join(", "), origin: "config" };
}

function inputSettings(command: AppCommand, config: SparseConfig): SparseConfig {
  return Object.fromEntries(
    commandSettings(command).settings.flatMap((spec) => {
      const value = config[spec.key];

      return spec.role === "input" && value !== undefined ? [[spec.key, value] as const] : [];
    }),
  );
}
