import type { AppCommand, CommandSettings } from "@neherlab/app-contracts";
import { settingCatalog } from "@neherlab/app-contracts/catalog";

class CatalogError extends Error {
  constructor(message: string) {
    super(message);
    this.name = "CatalogError";
  }
}

export const COMMANDS: readonly CommandSettings[] = settingCatalog.commands;

export function commandSettings(command: AppCommand): CommandSettings {
  const settings = COMMANDS.find((candidate) => candidate.command === command);

  if (settings === undefined) {
    throw new CatalogError(`the setting catalog has no command \`${command}\``);
  }

  return settings;
}
