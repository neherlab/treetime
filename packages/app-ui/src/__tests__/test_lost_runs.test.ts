import type { UiDraft } from "@neherlab/app-contracts";
import { ApiError } from "@neherlab/app-contracts/client";

import { forgetLostRuns } from "../analysis/lostRuns";
import { freshDraft } from "../store/draftSchema";
import { FakeServer, json, runRecord } from "./api_server";

const FRESH = freshDraft("clock");

const R1 = runRecord("r1", "created");

const R2 = runRecord("r2", "running");

const DRAFT: UiDraft = {
  ...FRESH,
  config: {
    ...FRESH.config,
    tree: "/treetime/runs/r1/inputs/tree.nwk",
    metadata: "/treetime/runs/r2/inputs/dates.tsv",
  },
  sources: {
    tree: { label: "tree.nwk", origin: "upload", run_id: "r1" },
    metadata: { label: "dates.tsv", origin: "run", run_id: "r2" },
  },
  upload_run_id: "r1",
};

describe("lost runs of a draft", () => {
  test("the files of a run that the back end no longer has leave the form", async () => {
    const server = new FakeServer({ "GET /api/runs/r2": () => json(R2) });

    const expected: UiDraft = {
      ...FRESH,
      config: { ...FRESH.config, metadata: "/treetime/runs/r2/inputs/dates.tsv" },
      sources: { metadata: { label: "dates.tsv", origin: "run", run_id: "r2" } },
    };

    await expect(forgetLostRuns(server.client(), DRAFT)).resolves.toStrictEqual({
      draft: expected,
      dropped: ["tree.nwk"],
      runs: new Map([["r2", R2]]),
    });
    expect(server.keys().toSorted()).toStrictEqual(["GET /api/runs/r1", "GET /api/runs/r2"]);
  });

  test("a draft whose runs all exist stays as it is", async () => {
    const server = new FakeServer({ "GET /api/runs/r1": () => json(R1), "GET /api/runs/r2": () => json(R2) });

    await expect(forgetLostRuns(server.client(), DRAFT)).resolves.toStrictEqual({
      draft: DRAFT,
      dropped: [],
      runs: new Map([
        ["r1", R1],
        ["r2", R2],
      ]),
    });
  });

  test("a draft without runs asks the back end nothing", async () => {
    const server = new FakeServer({});

    await expect(forgetLostRuns(server.client(), FRESH)).resolves.toStrictEqual({
      draft: FRESH,
      dropped: [],
      runs: new Map(),
    });
    expect(server.keys()).toStrictEqual([]);
  });

  test("a failed lookup is not taken for a lost run", async () => {
    const server = new FakeServer({
      "GET /api/runs/r1": () => json({ code: "internal_error", message: "the server is down", causes: [] }, 500),
      "GET /api/runs/r2": () => json(R2),
    });

    const error = await forgetLostRuns(server.client(), DRAFT).catch((failure: unknown) => failure);

    expect(error).toBeInstanceOf(ApiError);
    expect(error).toMatchObject({ status: 500, response: { code: "internal_error", message: "the server is down" } });
  });
});
