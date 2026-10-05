import type { AppCommand, SparseConfig, UiDraftSource } from "@neherlab/app-contracts";
import { configCheck } from "@neherlab/app-contracts/client";
import { useNavigate } from "@tanstack/react-router";
import { useCallback } from "react";

import { useApiContext } from "../api/context";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { normalizeConfig, settingValue } from "../settings/config";
import { baseName, pathList } from "../settings/inputs";
import { getAt, sameJson } from "../settings/json";
import { useDraftStore } from "../store/draft";

type ConfigLoadResult = { loaded: true; command: AppCommand } | { loaded: false; messages: string[] };

export function useConfigLoader() {
  const { client } = useApiContext();
  const navigate = useNavigate();

  return useCallback(
    async (
      text: string,
      fallbackCommand: AppCommand,
      keepInputs: boolean,
      folder?: string,
    ): Promise<ConfigLoadResult> => {
      const { draft, load } = useDraftStore.getState();
      const inputs = keepInputs ? inputSettings(draft.command, draft.config) : {};

      const { data: result } = await configCheck({
        client,
        body:
          folder === undefined
            ? { command: fallbackCommand, text, inputs }
            : { command: fallbackCommand, text, inputs, folder },
        throwOnError: true,
      });

      if (result.status === "invalid") {
        return { loaded: false, messages: result.messages };
      }

      const command = result.command;
      const specs = COMMAND_SETTINGS[command].settings;
      const config = normalizeConfig(specs, result.config);
      const sources: Record<string, UiDraftSource> = {};

      for (const spec of specs.filter((candidate) => candidate.role === "input")) {
        const value = settingValue(config, spec);
        const paths = pathList(value);
        const previous = draft.sources[spec.key];

        if (paths.length > 0) {
          sources[spec.key] =
            previous !== undefined && sameJson(value, getAt(draft.config, spec.path))
              ? previous
              : { label: paths.map(baseName).join(", "), origin: "config" };
        }
      }

      const { from_run_id: _fromRun, ...current } = draft;

      load({ ...current, command, config, sources });
      await navigate({ to: "/new" });

      return { loaded: true, command };
    },
    [client, navigate],
  );
}

function inputSettings(command: AppCommand, config: SparseConfig): SparseConfig {
  return Object.fromEntries(
    COMMAND_SETTINGS[command].settings.flatMap((spec) => {
      const value = getAt(config, spec.path);

      return spec.role === "input" && value !== undefined ? [[spec.key, value] as const] : [];
    }),
  );
}
