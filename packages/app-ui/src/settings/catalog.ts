import {
  openApiDocument,
  zSettingCatalog,
  type AppCommand,
  type CommandSettings,
  type SettingCatalog,
  type SettingSpec,
} from "@neherlab/app-contracts";

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
  const catalog: unknown = openApiDocument[SETTING_CATALOG_KEY];

  assertSettingCatalog(catalog);

  const entries = catalog.commands.map((command): [AppCommand, CommandSettings] => [command.command, command]);

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

function assertSettingCatalog(value: unknown): asserts value is SettingCatalog {
  zSettingCatalog.parse(value);
}
