import type { JobId } from "@neherlab/app-contracts";
import { ApiError, runsCancel, runsDelete, type ApiClient } from "@neherlab/app-contracts/client";

import { findRun } from "../analysis/pendingUpload";
import { runEventStream } from "../api/events";

export async function removeRun(client: ApiClient, id: JobId): Promise<void> {
  const record = await findRun(client, id);

  if (record === undefined) {
    return;
  }

  if (record.status === "running") {
    await runsCancel({ client, path: { id }, throwOnError: true });

    for await (const _event of runEventStream(client, id)) {
      continue;
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
