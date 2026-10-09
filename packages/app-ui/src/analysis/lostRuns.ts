import type { UiDraft } from "@neherlab/app-contracts";
import type { ApiClient } from "@neherlab/app-contracts/client";

import { forgetRuns, referencedRuns, type ForgottenRuns } from "../store/draftSchema";
import { findRun } from "./pendingUpload";

export async function forgetLostRuns(client: ApiClient, draft: UiDraft): Promise<ForgottenRuns> {
  const ids = referencedRuns(draft);
  const records = await Promise.all(ids.map(async (id) => findRun(client, id)));
  const lost = new Set(ids.filter((_, index) => records[index] === undefined));

  return lost.size === 0 ? { draft, dropped: [] } : forgetRuns(draft, lost);
}

export function lostFilesMessage(labels: readonly string[]): { title: string; description: string } {
  return {
    title: "Some input files are no longer available",
    description: `The back end no longer has ${labels.join(", ")}, so the form left them out. Add them again.`,
  };
}
