import type { AppCommand } from "@neherlab/app-contracts";
import { useNavigate } from "@tanstack/react-router";
import { useCallback } from "react";

import { useBridge } from "../BridgeContext";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { configCheckMessages } from "../settings/checks";
import { normalizeConfig, settingValue } from "../settings/config";
import { configTextCommand, withMissingInputs } from "../settings/configText";
import { baseName, pathList } from "../settings/inputs";
import { getAt, sameJson, zJsonObject } from "../settings/json";
import { useDraftStore, type InputSource } from "../store/draft";

export type ConfigLoadResult = { loaded: true; command: AppCommand } | { loaded: false; messages: string[] };

export function useConfigLoader() {
  const bridge = useBridge();
  const navigate = useNavigate();

  return useCallback(
    async (
      text: string,
      fallbackCommand: AppCommand,
      title: string | null,
      keepInputs: boolean,
    ): Promise<ConfigLoadResult> => {
      const command = configTextCommand(text) ?? fallbackCommand;
      const draft = useDraftStore.getState();
      const checked = keepInputs ? withMissingInputs(text, COMMAND_SETTINGS[command].specs, draft.config) : text;
      const result = await bridge.checkConfig({ command, text: checked });

      if (result.status === "invalid") {
        return { loaded: false, messages: configCheckMessages(result, false) };
      }

      const specs = COMMAND_SETTINGS[command].specs;
      const config = normalizeConfig(specs, zJsonObject.parse(result.config));
      const sources: Record<string, InputSource> = {};

      for (const spec of specs) {
        const value = settingValue(config, spec);

        const paths = pathList(value);
        const previous = draft.sources[spec.key];

        if (spec.pathRole === "input" && paths.length > 0) {
          sources[spec.key] =
            previous !== undefined && sameJson(value, getAt(draft.config, spec.path))
              ? previous
              : { label: paths.map(baseName).join(", "), origin: "config", size: null };
        }
      }

      useDraftStore.getState().load({
        command,
        config,
        sources,
        fromRunId: null,
        title: title ?? draft.title,
      });
      await navigate({ to: "/new" });

      return { loaded: true, command };
    },
    [bridge, navigate],
  );
}
