import { BridgeError, type RunRecordResult, type TreeTimeBridge } from "@neherlab/app-contracts";

export async function pendingUploadRun(
  bridge: Pick<TreeTimeBridge, "getRun">,
  id: string | null,
): Promise<RunRecordResult | undefined> {
  if (id === null) {
    return undefined;
  }

  try {
    const record = await bridge.getRun(id);

    return record.status === "created" ? record : undefined;
  } catch (error: unknown) {
    if (error instanceof BridgeError && error.response.code === "not_found") {
      return undefined;
    }

    throw error;
  }
}
