import type { RunRecord, SparseConfig } from "@neherlab/app-contracts";

import { commandSettings } from "./catalog";
import { normalizeConfig } from "./config";
import { runInputAssignments } from "./inputs";

const SEED_KEY = "seed";

interface RerunDraft {
  config: SparseConfig;
  inputLabels: Record<string, string>;
}

export function rerunDraft(record: RunRecord, settings: SparseConfig): RerunDraft {
  const specs = commandSettings(record.command).settings;
  const changed = new Set(record.changed_settings);
  const kept = new Set(specs.flatMap((spec) => (spec.role === "input" || changed.has(spec.key) ? [spec.key] : [])));

  const picked = Object.fromEntries(Object.entries(settings).filter(([key]) => kept.has(key)));
  const seeded = specs.some((spec) => spec.key === SEED_KEY) && record.seed !== undefined;
  const config = normalizeConfig(specs, seeded ? { ...picked, [SEED_KEY]: record.seed } : picked);

  const inputLabels = Object.fromEntries(
    runInputAssignments(record.command, record.inputs).map((input) => [input.key, input.label]),
  );

  return { config, inputLabels };
}
