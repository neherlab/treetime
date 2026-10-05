import type { AppCommand, UiDraft } from "@neherlab/app-contracts";

import { COMMAND_SETTINGS } from "../settings/catalog";
import { defaultConfig } from "../settings/config";

export function freshDraft(command: AppCommand): UiDraft {
  return {
    command,
    config: defaultConfig(COMMAND_SETTINGS[command].settings),
    sources: {},
    view: "main",
    search: "",
    changed_only: false,
    code_format: "cli",
  };
}
