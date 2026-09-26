import { ApiError, createApiClient } from "@neherlab/app-contracts/client";
import { describe, expect, test } from "vitest";

import { createWebSaveActions } from "../save-web";

const BASE_URL = "http://treetime.test";

function saver(routes: Record<string, () => Response>) {
  const requested: string[] = [];
  const saved: Array<[string, Blob]> = [];

  const fetch = (input: RequestInfo | URL, init?: RequestInit): Promise<Response> => {
    const request = new Request(input, init);
    const url = new URL(request.url);
    const key = `${request.method} ${url.pathname}${url.search}`;
    requested.push(key);

    return Promise.resolve((routes[key] ?? (() => new Response(null, { status: 500 })))());
  };

  const actions = createWebSaveActions(createApiClient({ baseUrl: BASE_URL, fetch }), (blob, name) => {
    saved.push([name, blob]);
  });

  return { actions, requested, saved };
}

describe("save_web", () => {
  test("saving a run file downloads the file bytes under the given name", async () => {
    const { actions, requested, saved } = saver({
      "GET /api/runs/r1/file?path=out%2Fclock.nwk": () =>
        new Response("(A,B);", { headers: { "Content-Type": "application/octet-stream" } }),
    });

    await expect(actions.saveRunFile("r1", "out/clock.nwk", "clock.nwk")).resolves.toBe(true);

    expect(requested).toStrictEqual(["GET /api/runs/r1/file?path=out%2Fclock.nwk"]);
    expect(await Promise.all(saved.map(async ([name, blob]) => [name, await blob.text()]))).toStrictEqual([
      ["clock.nwk", "(A,B);"],
    ]);
  });

  test("saving a run archive downloads the archive under the given name", async () => {
    const { actions, requested, saved } = saver({
      "GET /api/runs/r1/archive": () => new Response("PK", { headers: { "Content-Type": "application/zip" } }),
    });

    await expect(actions.saveRunArchive("r1", "clock.zip")).resolves.toBe(true);

    expect(requested).toStrictEqual(["GET /api/runs/r1/archive"]);
    expect(saved.map(([name]) => name)).toStrictEqual(["clock.zip"]);
  });

  test("a failed download rejects with the ApiError and saves nothing", async () => {
    const response = { code: "not_found", message: "no run with id `r9`", causes: [] };

    const { actions, saved } = saver({
      "GET /api/runs/r9/archive": () =>
        new Response(JSON.stringify(response), { status: 404, headers: { "Content-Type": "application/json" } }),
    });

    const error = await actions.saveRunArchive("r9", "r9.zip").catch((failure: unknown) => failure);

    expect(error).toBeInstanceOf(ApiError);
    expect(error).toMatchObject({ status: 404, response });
    expect(saved).toStrictEqual([]);
  });
});
