import { describe, expect, test } from "vitest";

import { router } from "../router";

const RUN_TABS = ["results", "settings", "log"] as const;

describe("run routes", () => {
  test.each(RUN_TABS)("the %s tab of a run is a child of the run page route", (tab) => {
    expect(router.matchRoutes(`/runs/r1/${tab}`).map((match) => match.routeId)).toStrictEqual([
      "__root__",
      "/runs/$id",
      `/runs/$id/${tab}`,
    ]);
  });

  test("every tab of one run shares the run page match, so the page and its event stream stay mounted", () => {
    const pageMatches = RUN_TABS.map(
      (tab) => router.matchRoutes(`/runs/r1/${tab}`).find((match) => match.routeId === "/runs/$id")?.id,
    );

    expect(new Set(pageMatches)).toStrictEqual(new Set([pageMatches[0]]));
    expect(pageMatches[0]).toBeDefined();
  });

  test("a bare run path matches the index route that opens the results tab", () => {
    expect(router.matchRoutes("/runs/r1").at(-1)?.routeId).toBe("/runs/$id/");
  });
});
