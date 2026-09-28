import type { RunRecord } from "@neherlab/app-contracts";

import { COMMAND_SETTINGS } from "./catalog";
import { normalizeConfig } from "./config";
import { runInputAssignments } from "./inputs";
import { zJsonObject, type JsonObject } from "./json";

interface RerunDraft {
  config: JsonObject;
  inputLabels: Record<string, string>;
}

export function rerunDraft(record: RunRecord): RerunDraft {
  const specs = COMMAND_SETTINGS[record.command].specs;

  return {
    config: normalizeConfig(specs, zJsonObject.parse(record.config)),
    inputLabels: Object.fromEntries(
      runInputAssignments(record.command, record.inputs).map((input) => [input.key, input.label]),
    ),
  };
}
