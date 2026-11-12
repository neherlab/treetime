import { hasRunView, runViews } from "../views";

describe("views of a run", () => {
  test.each([
    ["timetree", ["Overview", "Clock", "Population", "Model", "Log", "Settings"]],
    ["clock", ["Overview", "Log", "Settings"]],
    ["ancestral", ["Overview", "Model", "Log", "Settings"]],
    ["optimize", ["Overview", "Model", "Log", "Settings"]],
    ["prune", ["Overview", "Log", "Settings"]],
  ] as const)("a finished %s run has the views %j", (command, expected) => {
    expect(runViews(command, "ok").map((view) => view.label)).toStrictEqual(expected);
  });

  test.each(["created", "running", "error", "cancelled", "interrupted"] as const)(
    "a %s time tree run has only the overview, the log, and the settings",
    (status) => {
      expect(runViews("timetree", status).map((view) => view.key)).toStrictEqual(["overview", "log", "settings"]);
    },
  );

  test("each view has its own route", () => {
    expect(runViews("timetree", "ok").map((view) => view.to)).toStrictEqual([
      "/runs/$id/results",
      "/runs/$id/clock",
      "/runs/$id/population",
      "/runs/$id/model",
      "/runs/$id/log",
      "/runs/$id/settings",
    ]);
  });

  test("a view that the analysis lacks is not available", () => {
    expect([hasRunView("clock", "ok", "population"), hasRunView("timetree", "ok", "population")]).toStrictEqual([
      false,
      true,
    ]);
  });
});
