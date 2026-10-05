import type { RunRecord } from "@neherlab/app-contracts";
import { ApiError, runsGet, type ApiClient } from "@neherlab/app-contracts/client";

export async function pendingUploadRun(client: ApiClient, id: string | undefined): Promise<RunRecord | undefined> {
  if (id === undefined) {
    return undefined;
  }

  try {
    const { data: record } = await runsGet({ client, path: { id }, throwOnError: true });

    return record.status === "created" ? record : undefined;
  } catch (error: unknown) {
    if (error instanceof ApiError && error.response.code === "not_found") {
      return undefined;
    }

    throw error;
  }
}
