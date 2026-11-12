import type { AppCommand } from "@neherlab/app-contracts";

import { freshDraft } from "../draftSchema";

export function blankDraftState(command: AppCommand) {
  return {
    draft: freshDraft(command),
    viewing: undefined,
    lastShown: undefined,
    undo: undefined,
    queuedExamples: [],
    selections: {},
    uploadRunCreation: undefined,
    busyRuns: [],
  };
}
