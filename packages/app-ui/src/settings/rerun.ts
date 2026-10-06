import type { AppCommand, RunRecord, SparseConfig } from "@neherlab/app-contracts";

import { commandSettings } from "./catalog";
import { carryOverDraft } from "./commands";
import { normalizeConfig } from "./config";
import { runInputAssignments } from "./inputs";

interface RerunDraft {
  config: SparseConfig;
  inputLabels: Record<string, string>;
}

export function rerunDraft(record: RunRecord, command: AppCommand = record.command): RerunDraft {
  const config = normalizeConfig(commandSettings(record.command).settings, record.config);

  const inputLabels = Object.fromEntries(
    runInputAssignments(record.command, record.inputs).map((input) => [input.key, input.label]),
  );

  if (command === record.command) {
    return { config, inputLabels };
  }

  const carried = carryOverDraft(command, record.command, config, inputLabels);

  return { config: carried.config, inputLabels: carried.sources };
}
