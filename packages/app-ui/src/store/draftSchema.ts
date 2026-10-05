import type { AppCommand, UiDraft } from "@neherlab/app-contracts";

import { commandSettings } from "../settings/catalog";
import { defaultConfig } from "../settings/config";

export function freshDraft(command: AppCommand): UiDraft {
  return {
    command,
    config: defaultConfig(commandSettings(command).settings),
    sources: {},
    view: "main",
    search: "",
    changed_only: false,
    code_format: "cli",
  };
}
