import {
  openApiDocument,
  zSettingCatalog,
  type AppCommand,
  type zCommandSettings,
  type zSettingSpec,
} from "@neherlab/app-contracts";
import type * as z from "zod";

import { zJsonValue, type JsonValue } from "./json";

export type SettingSpec = Omit<z.infer<typeof zSettingSpec>, "default_value" | "examples"> & {
  default_value: JsonValue;
  examples: JsonValue[];
};

export type CommandSettings = Omit<z.infer<typeof zCommandSettings>, "settings"> & { specs: SettingSpec[] };

const SETTING_CATALOG_KEY = "x-setting-catalog";

class CatalogError extends Error {
  constructor(message: string) {
    super(message);
    this.name = "CatalogError";
  }
}

export const COMMAND_SETTINGS: Readonly<Record<AppCommand, CommandSettings>> = commandSettingsByCommand();

export function groupedSpecs(settings: CommandSettings, specs: readonly SettingSpec[]): Array<[string, SettingSpec[]]> {
  return settings.groups.flatMap((group) => {
    const members = specs.filter((spec) => spec.group === group);

    return members.length > 0 ? [[group, members]] : [];
  });
}

function commandSettingsByCommand(): Record<AppCommand, CommandSettings> {
  const catalog = zSettingCatalog.parse(openApiDocument[SETTING_CATALOG_KEY]);

  const entries = catalog.commands.map((command): [AppCommand, CommandSettings] => [
    command.command,
    {
      command: command.command,
      inputs: command.inputs,
      uses_dates: command.uses_dates,
      groups: command.groups,
      specs: command.settings.map((spec) => ({
        ...spec,
        default_value: zJsonValue.parse(spec.default_value),
        examples: spec.examples.map((example) => zJsonValue.parse(example)),
      })),
    },
  ]);

  const byCommand = new Map(entries);

  const lookup = (command: AppCommand): CommandSettings => {
    const settings = byCommand.get(command);

    if (settings === undefined) {
      throw new CatalogError(`the setting catalog has no command \`${command}\``);
    }

    return settings;
  };

  return {
    timetree: lookup("timetree"),
    clock: lookup("clock"),
    ancestral: lookup("ancestral"),
    mugration: lookup("mugration"),
    optimize: lookup("optimize"),
    prune: lookup("prune"),
  };
}
