import type { RunRecordResult } from "@neherlab/app-contracts";

import { COMMAND_SETTINGS } from "./catalog";
import { normalizeConfig } from "./config";
import { runInputAssignments } from "./inputs";
import { zJsonObject, type JsonObject } from "./json";

export interface RerunDraft {
  config: JsonObject;
  inputLabels: Record<string, string>;
  title: string;
}

export function rerunDraft(record: RunRecordResult): RerunDraft {
  const specs = COMMAND_SETTINGS[record.command].specs;

  return {
    config: normalizeConfig(specs, zJsonObject.parse(record.config)),
    inputLabels: Object.fromEntries(
      runInputAssignments(record.command, record.inputs).map((input) => [input.key, input.label]),
    ),
    title: `${record.title} (edited)`,
  };
}
