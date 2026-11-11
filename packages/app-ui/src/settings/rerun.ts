import type { RunRecord, SparseConfig } from "@neherlab/app-contracts";

import { commandSettings } from "./catalog";
import { normalizeConfig } from "./config";
import { runInputAssignments } from "./inputs";

interface RerunDraft {
  config: SparseConfig;
  inputLabels: Record<string, string>;
}

export function rerunDraft(record: RunRecord): RerunDraft {
  const specs = commandSettings(record.command).settings;
  const changed = new Set(record.changed_settings);
  const kept = new Set(specs.flatMap((spec) => (spec.role === "input" || changed.has(spec.key) ? [spec.key] : [])));
  const config = normalizeConfig(
    specs,
    Object.fromEntries(Object.entries(record.config).filter(([key]) => kept.has(key))),
  );

  const inputLabels = Object.fromEntries(
    runInputAssignments(record.command, record.inputs).map((input) => [input.key, input.label]),
  );

  return { config, inputLabels };
}
