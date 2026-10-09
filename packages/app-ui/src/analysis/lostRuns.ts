import type { JobId, RunRecord, UiDraft } from "@neherlab/app-contracts";
import type { ApiClient } from "@neherlab/app-contracts/client";

import { forgetRuns, referencedRuns, type ForgottenRuns } from "../store/draftSchema";
import { findRun } from "./pendingUpload";

export async function forgetLostRuns(client: ApiClient, draft: UiDraft): Promise<CheckedRuns> {
  const ids = referencedRuns(draft);
  const records = await Promise.all(ids.map(async (id) => findRun(client, id)));
  const runs = new Map(records.flatMap((record) => (record === undefined ? [] : [[record.id, record] as const])));
  const lost = new Set(ids.filter((id) => !runs.has(id)));

  return { ...(lost.size === 0 ? { draft, dropped: [] } : forgetRuns(draft, lost)), runs };
}

export function lostFilesMessage(labels: readonly string[]) {
  return {
    title: "Some input files are no longer available",
    description: `The back end no longer has ${labels.join(", ")}, so the form left them out. Add them again.`,
  };
}

export interface CheckedRuns extends ForgottenRuns {
  runs: ReadonlyMap<JobId, RunRecord>;
}
