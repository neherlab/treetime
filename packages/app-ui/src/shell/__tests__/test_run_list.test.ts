import type { RunSummary } from "@neherlab/app-contracts";

import { groupedRunOrder, groupRuns, listedRuns, type RunGroup } from "../runList";

function run(id: string, overrides: Partial<RunSummary>): RunSummary {
  return {
    id,
    title: id,
    command: "timetree",
    status: "ok",
    pinned: false,
    group_id: id,
    size: 0,
    created_at: "2026-09-25T08:00:00Z",
    changed_settings: [],
    headline: {},
    package_name: `${id}.zip`,
    log_name: `${id}.log.txt`,
    config_name: `${id}.yaml`,
    ...overrides,
  };
}

describe("run list", () => {
  test("drafts waiting for uploads are not listed", () => {
    const runs = [run("a", { status: "created" }), run("b", {})];

    expect(listedRuns(runs, "", null).map((listed) => listed.id)).toStrictEqual(["b"]);
  });

  test("the filter matches every word over title and analysis", () => {
    const runs = [
      run("a", { title: "ebola/100, time tree" }),
      run("b", { title: "ebola/100, time tree: skyline" }),
      run("c", { title: "ebola/100, clock", command: "clock" }),
    ];

    expect(listedRuns(runs, "time skyline", null).map((listed) => listed.id)).toStrictEqual(["b"]);
  });

  test("hidden runs are not listed", () => {
    const runs = [run("a", {}), run("b", {})];

    expect(listedRuns(runs, "", null, new Set(["a"])).map((listed) => listed.id)).toStrictEqual(["b"]);
  });

  test("the command filter keeps one command", () => {
    const runs = [run("a", {}), run("b", { command: "clock" })];

    expect(listedRuns(runs, "", "clock").map((listed) => listed.id)).toStrictEqual(["b"]);
  });

  test("each derivation group lists its variants under their source, newest group first", () => {
    const runs = [
      run("v2", { group_id: "a", source_run_id: "a", created_at: "2026-09-25T10:00:00Z" }),
      run("b", { group_id: "b", created_at: "2026-09-25T09:00:00Z" }),
      run("v1", { group_id: "a", source_run_id: "a", created_at: "2026-09-25T08:30:00Z" }),
      run("v11", { group_id: "a", source_run_id: "v1", created_at: "2026-09-25T09:30:00Z" }),
      run("a", { group_id: "a", created_at: "2026-09-25T08:00:00Z" }),
    ];

    expect(groupedIds(groupRuns(runs))).toStrictEqual([
      ["a", ["a:0", "v1:1", "v11:2", "v2:1"]],
      ["b", ["b:0"]],
    ]);
  });

  test("groups with a pinned run come first", () => {
    const runs = [
      run("a", { created_at: "2026-09-25T10:00:00Z" }),
      run("b", { pinned: true, created_at: "2026-09-20T08:00:00Z" }),
    ];

    expect(groupedIds(groupRuns(runs))).toStrictEqual([
      ["b", ["b:0"]],
      ["a", ["a:0"]],
    ]);
  });

  test("a variant whose source is not listed starts its own line in the group", () => {
    const runs = [
      run("v1", { group_id: "a", source_run_id: "a", created_at: "2026-09-25T09:00:00Z" }),
      run("v2", { group_id: "a", source_run_id: "x", created_at: "2026-09-25T10:00:00Z" }),
    ];

    expect(groupedIds(groupRuns(runs))).toStrictEqual([["a", ["v1:0", "v2:0"]]]);
  });

  test("the keyboard order follows the groups", () => {
    const runs = [
      run("b", { created_at: "2026-09-25T09:00:00Z" }),
      run("v", { group_id: "a", source_run_id: "a", created_at: "2026-09-25T10:00:00Z" }),
      run("a", { created_at: "2026-09-25T08:00:00Z" }),
    ];

    expect(groupedRunOrder(groupRuns(runs)).map((ordered) => ordered.id)).toStrictEqual(["a", "v", "b"]);
  });
});

function groupedIds(groups: readonly RunGroup[]): Array<[string, string[]]> {
  return groups.map((group) => [group.id, group.entries.map((entry) => `${entry.run.id}:${entry.depth}`)]);
}
