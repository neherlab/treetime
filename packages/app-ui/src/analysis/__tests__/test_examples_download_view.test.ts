import type { DataIndex, ExamplesDownloadStatus } from "@neherlab/app-contracts";

import type { ExampleIndex } from "../exampleData";
import { examplesDownloadView } from "../examplesDownloadView";

const INDEX: DataIndex = {
  files: [{ path: "zika/20/tree.nwk", size: 11 }],
  datasets: [{ name: "zika/20", inputs: [{ kind: "tree", path: "zika/20/tree.nwk" }] }],
  configs: [],
  pipelines: [],
};

const EMPTY: ExampleIndex = { origin: "local", folder: "/app/examples", index: undefined };

const FILLED: ExampleIndex = { origin: "local", folder: "/app/examples", index: INDEX };

const WEB: ExampleIndex = { origin: "web", index: INDEX };

describe("examples download view", () => {
  test("a local app whose examples folder has no index offers the download", () => {
    expect(examplesDownloadView(EMPTY, { download: { state: "idle" } })).toStrictEqual({ kind: "offer" });
  });

  test("the web apps and an examples folder with an index show nothing", () => {
    expect([examplesDownloadView(WEB, undefined), examplesDownloadView(FILLED, undefined)]).toStrictEqual([
      { kind: "hidden" },
      { kind: "hidden" },
    ]);
  });

  test("a running download shows its share of the archive", () => {
    const status: ExamplesDownloadStatus = { download: { state: "running", received: 250, total: 1000 }, seq: 3 };

    expect(examplesDownloadView(EMPTY, status)).toStrictEqual({ kind: "running", received: 250, percent: 25 });
  });

  test("a running download of unknown size shows no share", () => {
    const status: ExamplesDownloadStatus = { download: { state: "running", received: 250 }, seq: 3 };

    expect(examplesDownloadView(EMPTY, status)).toStrictEqual({
      kind: "running",
      received: 250,
      percent: undefined,
    });
  });

  test("a failed download shows its message and offers a retry", () => {
    const status: ExamplesDownloadStatus = { download: { state: "failed", message: "HTTP 404" }, seq: 4 };

    expect(examplesDownloadView(EMPTY, status)).toStrictEqual({ kind: "failed", message: "HTTP 404" });
  });

  test("a finished download hides the offer once the folder has an index", () => {
    const status: ExamplesDownloadStatus = { download: { state: "done", received: 1000 }, seq: 5 };

    expect(examplesDownloadView(FILLED, status)).toStrictEqual({ kind: "hidden" });
  });
});
