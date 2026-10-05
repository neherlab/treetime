import { ApiError } from "@neherlab/app-contracts/client";
import { describe, expect, test } from "vitest";

import { pendingUploadRun } from "../analysis/pendingUpload";
import { FakeServer, json, RECORD } from "./api_server";

const { started_at: _startedAt, ...RUNNING } = RECORD;

const CREATED = { ...RUNNING, title: "Uploaded inputs", config: {}, status: "created" };

describe("pending upload run", () => {
  test("a created run is the pending upload run", async () => {
    const server = new FakeServer({ "GET /api/runs/r1": () => json(CREATED) });

    await expect(pendingUploadRun(server.client(), "r1")).resolves.toStrictEqual(CREATED);
  });

  test("a started run or a missing run is no pending upload run", async () => {
    const server = new FakeServer({ "GET /api/runs/r2": () => json({ ...CREATED, id: "r2", status: "running" }) });
    const client = server.client();

    await expect(
      Promise.all([
        pendingUploadRun(client, "r2"),
        pendingUploadRun(client, "r9"),
        pendingUploadRun(client, undefined),
      ]),
    ).resolves.toStrictEqual([undefined, undefined, undefined]);
    expect(server.keys().toSorted()).toStrictEqual(["GET /api/runs/r2", "GET /api/runs/r9"]);
  });

  test("a failed lookup is not taken for a missing run", async () => {
    const server = new FakeServer({
      "GET /api/runs/r1": () => json({ code: "internal_error", message: "the server is down", causes: [] }, 500),
    });

    const error = await pendingUploadRun(server.client(), "r1").catch((failure: unknown) => failure);

    expect(error).toBeInstanceOf(ApiError);
    expect(error).toMatchObject({ status: 500, response: { code: "internal_error", message: "the server is down" } });
  });
});
