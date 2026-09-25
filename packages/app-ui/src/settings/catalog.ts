import type { AppCommand } from "@neherlab/app-contracts";

import { commandSettings, type CommandSettings } from "./schema";

export const COMMAND_SETTINGS: Readonly<Record<AppCommand, CommandSettings>> = {
  timetree: commandSettings("timetree"),
  clock: commandSettings("clock"),
  ancestral: commandSettings("ancestral"),
  mugration: commandSettings("mugration"),
  optimize: commandSettings("optimize"),
  prune: commandSettings("prune"),
};
