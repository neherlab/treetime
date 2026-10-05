import type { ExamplesDownloadStatus } from "@neherlab/app-contracts";
import { describe, expect, test } from "vitest";

import { newerDownloadStatus } from "../examplesDownload";

const RUNNING_AT_5: ExamplesDownloadStatus = { download: { state: "running", received: 10, total: 100 }, seq: 5 };

const DONE_AT_7: ExamplesDownloadStatus = { download: { state: "done", received: 100 }, seq: 7 };

const IDLE: ExamplesDownloadStatus = { download: { state: "idle" } };

describe("examples download status order", () => {
  test("an event newer than the snapshot replaces it", () => {
    expect(newerDownloadStatus(RUNNING_AT_5, DONE_AT_7)).toBe(DONE_AT_7);
  });

  test("a slow snapshot older than the received event does not overwrite it", () => {
    expect(newerDownloadStatus(DONE_AT_7, RUNNING_AT_5)).toBe(DONE_AT_7);
  });

  test("a status without an event loses to any reported event", () => {
    expect([newerDownloadStatus(RUNNING_AT_5, IDLE), newerDownloadStatus(IDLE, RUNNING_AT_5)]).toStrictEqual([
      RUNNING_AT_5,
      RUNNING_AT_5,
    ]);
  });

  test("a missing side keeps the other", () => {
    expect([newerDownloadStatus(undefined, IDLE), newerDownloadStatus(DONE_AT_7, undefined)]).toStrictEqual([
      IDLE,
      DONE_AT_7,
    ]);
  });
});
