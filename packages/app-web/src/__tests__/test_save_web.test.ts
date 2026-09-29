import { createApiClient } from "@neherlab/app-contracts/client";
import { describe, expect, test } from "vitest";

import { createWebSaveActions } from "../save-web";

const BASE_URL = "http://treetime.test";

function saver() {
  const requested: string[] = [];
  const saved: Array<[string, string]> = [];

  const fetch = (input: RequestInfo | URL): Promise<Response> => {
    requested.push(new Request(input).url);

    return Promise.resolve(new Response(null, { status: 500 }));
  };

  const actions = createWebSaveActions(createApiClient({ baseUrl: BASE_URL, fetch }), (url, name) => {
    saved.push([url, name]);
  });

  return { actions, requested, saved };
}

describe("save_web", () => {
  test("saving a run file hands the file URL to the browser under the given name", async () => {
    const { actions, requested, saved } = saver();

    await expect(actions.saveRunFile("r1", "out/clock.nwk", "clock.nwk")).resolves.toBe(true);

    expect({ requested, saved }).toStrictEqual({
      requested: [],
      saved: [["http://treetime.test/api/runs/r1/file?path=out%2Fclock.nwk", "clock.nwk"]],
    });
  });

  test("saving a run archive hands the archive URL to the browser under the given name", async () => {
    const { actions, requested, saved } = saver();

    await expect(actions.saveRunArchive("r1", "clock.zip")).resolves.toBe(true);

    expect({ requested, saved }).toStrictEqual({
      requested: [],
      saved: [["http://treetime.test/api/runs/r1/archive", "clock.zip"]],
    });
  });
});
