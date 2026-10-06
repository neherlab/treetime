import type { AppCommand, RunRecord, SparseConfig } from "@neherlab/app-contracts";

import { commandSettings } from "./catalog";
import { carryOverConfig } from "./config";

type ConfigOf<C extends AppCommand> = Extract<RunRecord, { command: C }>["config"];

export type SettingKey<C extends AppCommand> = keyof ConfigOf<C> & string;

export const COMMAND_VERBS: Record<AppCommand, string> = {
  timetree: "Run time tree",
  clock: "Run clock check",
  ancestral: "Reconstruct sequences",
  homoplasy: "Find homoplasies",
  mugration: "Reconstruct traits",
  optimize: "Optimize branch lengths",
  prune: "Prune tree",
};

export function carryOverDraft<Source>(
  target: AppCommand,
  source: AppCommand,
  config: SparseConfig,
  sources: Readonly<Record<string, Source>>,
): CarriedDraft<Source> {
  const targetSpecs = commandSettings(target).settings;
  const keys = new Set(targetSpecs.map((spec) => spec.key));

  return {
    config: carryOverConfig(targetSpecs, commandSettings(source).settings, config),
    sources: Object.fromEntries(Object.entries(sources).filter(([key]) => keys.has(key))),
  };
}

export function commandSwitchNote(requested: AppCommand, loaded: AppCommand): string | undefined {
  return requested === loaded ? undefined : `Switched to ${commandSettings(loaded).title}`;
}

export interface CarriedDraft<Source> {
  config: SparseConfig;
  sources: Record<string, Source>;
}
