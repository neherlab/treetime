import type { JobId } from "@neherlab/app-contracts";

import { useDraftStore } from "../store/draft";

export async function whileRunBusy<T>(runId: JobId, work: () => Promise<T>): Promise<T> {
  useDraftStore.getState().beginRunWork(runId);

  try {
    return await work();
  } finally {
    useDraftStore.getState().endRunWork(runId);
  }
}
