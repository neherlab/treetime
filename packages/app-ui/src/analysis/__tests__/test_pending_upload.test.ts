import { BridgeError, type RunRecordResult } from "@neherlab/app-contracts";
import { describe, expect, test } from "vitest";

import { pendingUploadRun } from "../pendingUpload";

const CREATED: RunRecordResult = {
  id: "r1",
  title: "Uploaded inputs",
  command: "clock",
  config: {},
  status: "created",
  pinned: false,
  created_at: "2026-09-25T10:00:00Z",
  started_at: null,
  finished_at: null,
  duration_seconds: null,
  treetime_version: "1.0.0",
  inputs: [],
  config_hash: null,
  changed_settings: [],
  headline: {},
  output_files: [],
  error: null,
};

describe("pending upload run", () => {
  test("a created run is the pending upload run", async () => {
    await expect(pendingUploadRun({ getRun: () => Promise.resolve(CREATED) }, "r1")).resolves.toStrictEqual(CREATED);
  });

  test("a started run or a missing run is no pending upload run", async () => {
    const missing = new BridgeError({ code: "not_found", message: "no run with id `r1`", causes: [] });

    await expect(
      Promise.all([
        pendingUploadRun({ getRun: () => Promise.resolve({ ...CREATED, status: "running" }) }, "r1"),
        pendingUploadRun({ getRun: () => Promise.reject(missing) }, "r1"),
        pendingUploadRun({ getRun: () => Promise.resolve(CREATED) }, null),
      ]),
    ).resolves.toStrictEqual([undefined, undefined, undefined]);
  });

  test("a failed lookup is not taken for a missing run", async () => {
    const failure = new BridgeError({ code: "internal_error", message: "the server is down", causes: [] });

    await expect(pendingUploadRun({ getRun: () => Promise.reject(failure) }, "r1")).rejects.toBe(failure);
  });
});
