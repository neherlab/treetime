import type { AppCommand, RunRecord } from "@neherlab/app-contracts";

import { commandSettings } from "./catalog";

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

export function commandSwitchNote(requested: AppCommand, loaded: AppCommand): string | undefined {
  return requested === loaded ? undefined : `Switched to ${commandSettings(loaded).title}`;
}
