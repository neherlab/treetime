import type { RunSummary } from "@neherlab/app-contracts";

import { duplicateRun } from "../duplicateRun";

const RUN: RunSummary = {
  id: "a",
  title: "zika-20-time-tree",
  command: "timetree",
  status: "ok",
  pinned: false,
  group_id: "a",
  size: 0,
  created_at: "2026-11-12T08:00:00Z",
  config_hash: "h",
  seed: 4_242,
  changed_settings: [],
  headline: {},
  package_name: "zika-20-time-tree-a.zip",
};

describe("duplicate run of a draft", () => {
  test("a draft without a seed duplicates a run with the same hash, and a new run draws another seed", () => {
    expect(duplicateRun([RUN], "timetree", { config_hash: "h" })).toStrictEqual({ run: RUN, repeatsSeed: false });
  });

  test("a draft with the seed of the run duplicates it and repeats its random step", () => {
    expect(duplicateRun([RUN], "timetree", { config_hash: "h", seed: 4_242 })).toStrictEqual({
      run: RUN,
      repeatsSeed: true,
    });
  });

  test("a draft with another seed does not duplicate the run", () => {
    expect(duplicateRun([RUN], "timetree", { config_hash: "h", seed: 7 })).toBeUndefined();
  });

  test("a draft with another hash, another command, or an unfinished run has no duplicate", () => {
    expect([
      duplicateRun([RUN], "timetree", { config_hash: "other" }),
      duplicateRun([RUN], "clock", { config_hash: "h" }),
      duplicateRun([{ ...RUN, status: "error" }], "timetree", { config_hash: "h" }),
      duplicateRun([RUN], "timetree", {}),
      duplicateRun([RUN], "timetree", undefined),
    ]).toStrictEqual([undefined, undefined, undefined, undefined, undefined]);
  });
});
