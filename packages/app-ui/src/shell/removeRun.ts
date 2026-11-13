import type { JobId } from "@neherlab/app-contracts";
import { ApiError, runsCancel, runsDelete, type ApiClient } from "@neherlab/app-contracts/client";

import { findRun } from "../analysis/pendingUpload";
import { runEventStream } from "../api/events";

export async function removeRun(client: ApiClient, id: JobId, deadline: AbortSignal): Promise<void> {
  const record = await findRun(client, id);

  if (record === undefined) {
    return;
  }

  if (record.status === "running") {
    const { data: cancel } = await runsCancel({ client, path: { id }, throwOnError: true });

    try {
      for await (const _event of runEventStream(client, id, { from: cancel.from ?? 0, signal: deadline })) {
        continue;
      }
    } catch (error: unknown) {
      if (!deadline.aborted) {
        throw error;
      }
    }

    if (deadline.aborted) {
      throw new Error(`${record.title} did not stop after the cancel request; delete it again once it has ended`);
    }
  }

  try {
    await runsDelete({ client, path: { id }, throwOnError: true });
  } catch (error: unknown) {
    if (!(error instanceof ApiError && error.response.code === "not_found")) {
      throw error;
    }
  }
}
