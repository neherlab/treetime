import { router } from "../router";

const RUN_TABS = ["results", "settings", "log"] as const;

describe("root route", () => {
  test("the root path matches the index route that opens a new analysis", () => {
    expect(router.matchRoutes("/").at(-1)?.routeId).toBe("/");
  });
});

describe("run routes", () => {
  test.each([...RUN_TABS])("the %s tab of a run is a child of the run page route in the workspace", (tab) => {
    expect(router.matchRoutes(`/runs/r1/${tab}`).map((match) => match.routeId)).toStrictEqual([
      "__root__",
      "/workspace",
      "/workspace/runs/$id",
      `/workspace/runs/$id/${tab}`,
    ]);
  });

  test("every tab of one run shares the run page match, so the page and its event stream stay mounted", () => {
    const pageMatches = RUN_TABS.map(
      (tab) => router.matchRoutes(`/runs/r1/${tab}`).find((match) => match.routeId === "/workspace/runs/$id")?.id,
    );

    expect(new Set(pageMatches)).toStrictEqual(new Set([pageMatches[0]]));
    expect(pageMatches[0]).toBeDefined();
  });

  test("a bare run path matches the index route that opens the results tab", () => {
    expect(router.matchRoutes("/runs/r1").at(-1)?.routeId).toBe("/workspace/runs/$id/");
  });
});

describe("workspace route", () => {
  test("a new analysis and the runs share one workspace match, so the options stay mounted between them", () => {
    const workspaceMatches = ["/new", "/runs/r1/results", "/runs/r2/log"].map(
      (path) => router.matchRoutes(path).find((match) => match.routeId === "/workspace")?.id,
    );

    expect(new Set(workspaceMatches)).toStrictEqual(new Set([workspaceMatches[0]]));
    expect(workspaceMatches[0]).toBeDefined();
  });

  test("the comparison of two runs stays outside the workspace", () => {
    expect(router.matchRoutes("/compare/r1/r2").map((match) => match.routeId)).toStrictEqual([
      "__root__",
      "/compare/$a/$b",
    ]);
  });
});
