import type { RunSummary } from "@neherlab/app-contracts";
import { DateTime } from "luxon";

import { changedFlags, groupRuns, listedRuns } from "../runList";

function run(id: string, overrides: Partial<RunSummary>): RunSummary {
  return {
    id,
    title: id,
    command: "timetree",
    status: "ok",
    pinned: false,
    created_at: "2026-09-25T08:00:00Z",
    changed_settings: [],
    headline: {},
    ...overrides,
  };
}

const NOW = DateTime.fromISO("2026-09-25T12:00:00Z", { zone: "utc" });

describe("run list", () => {
  test("drafts waiting for uploads are not listed", () => {
    const runs = [run("a", { status: "created" }), run("b", {})];

    expect(listedRuns(runs, "", null).map((listed) => listed.id)).toStrictEqual(["b"]);
  });

  test("the filter matches every word over title, command and changed flags", () => {
    const runs = [
      run("a", { title: "Baseline" }),
      run("b", { title: "Skyline", changed_settings: ["coalescent"] }),
      run("c", { title: "Clock check", command: "clock" }),
    ];

    expect(listedRuns(runs, "time --coalescent", null).map((listed) => listed.id)).toStrictEqual(["b"]);
  });

  test("the command filter keeps one command", () => {
    const runs = [run("a", {}), run("b", { command: "clock" })];

    expect(listedRuns(runs, "", "clock").map((listed) => listed.id)).toStrictEqual(["b"]);
  });

  test("pinned runs come first, then one group per day", () => {
    const runs = [
      run("a", {}),
      run("b", { pinned: true, created_at: "2026-09-20T08:00:00Z" }),
      run("c", { created_at: "2026-09-24T08:00:00Z" }),
    ];

    expect(groupRuns(runs, NOW).map((group) => [group.label, group.runs.map((grouped) => grouped.id)])).toStrictEqual([
      ["Pinned", ["b"]],
      ["Today", ["a"]],
      ["Yesterday", ["c"]],
    ]);
  });

  test("changed settings show as flags", () => {
    expect(changedFlags(run("a", { command: "clock", changed_settings: ["keep_root", "clock_filter"] }))).toStrictEqual(
      ["--keep-root", "--clock-filter"],
    );
  });
});
