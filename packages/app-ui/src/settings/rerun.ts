import type { RunRecord, SparseConfig } from "@neherlab/app-contracts";

import { COMMAND_SETTINGS } from "./catalog";
import { normalizeConfig } from "./config";
import { runInputAssignments } from "./inputs";

interface RerunDraft {
  config: SparseConfig;
  inputLabels: Record<string, string>;
}

export function rerunDraft(record: RunRecord): RerunDraft {
  const specs = COMMAND_SETTINGS[record.command].settings;

  return {
    config: normalizeConfig(specs, record.config),
    inputLabels: Object.fromEntries(
      runInputAssignments(record.command, record.inputs).map((input) => [input.key, input.label]),
    ),
  };
}
