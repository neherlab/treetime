import type { DatasetCatalog, ExamplesDownloadStatus } from "@neherlab/app-contracts";
import { describe, expect, test } from "vitest";

import { examplesDownloadView } from "../examplesDownloadView";

const EMPTY: DatasetCatalog = { datasets: [], examples: [] };

const FILLED: DatasetCatalog = {
  datasets: [{ name: "zika/20", files: ["tree.nwk"], inputs: [] }],
  examples: [],
};

describe("examples download view", () => {
  test("a local app with an empty examples folder offers the download", () => {
    expect(examplesDownloadView(true, EMPTY, { download: { state: "idle" } })).toStrictEqual({ kind: "offer" });
  });

  test("the web app and a filled examples folder show nothing", () => {
    expect([
      examplesDownloadView(false, EMPTY, undefined),
      examplesDownloadView(true, FILLED, undefined),
    ]).toStrictEqual([{ kind: "hidden" }, { kind: "hidden" }]);
  });

  test("a running download shows its share of the archive", () => {
    const status: ExamplesDownloadStatus = { download: { state: "running", received: 250, total: 1000 }, seq: 3 };

    expect(examplesDownloadView(true, EMPTY, status)).toStrictEqual({ kind: "running", received: 250, percent: 25 });
  });

  test("a running download of unknown size shows no share", () => {
    const status: ExamplesDownloadStatus = { download: { state: "running", received: 250 }, seq: 3 };

    expect(examplesDownloadView(true, EMPTY, status)).toStrictEqual({
      kind: "running",
      received: 250,
      percent: undefined,
    });
  });

  test("a failed download shows its message and offers a retry", () => {
    const status: ExamplesDownloadStatus = { download: { state: "failed", message: "HTTP 404" }, seq: 4 };

    expect(examplesDownloadView(true, EMPTY, status)).toStrictEqual({ kind: "failed", message: "HTTP 404" });
  });

  test("a finished download hides the offer once the catalog lists the datasets", () => {
    const status: ExamplesDownloadStatus = { download: { state: "done", received: 1000 }, seq: 5 };

    expect(examplesDownloadView(true, FILLED, status)).toStrictEqual({ kind: "hidden" });
  });
});
