import type { JobId, RunStatus } from "@neherlab/app-contracts";

import { comparePair } from "../comparePair";

const RUNS: Array<{ id: JobId; status: RunStatus }> = [
  { id: "r1", status: "ok" },
  { id: "r2", status: "ok" },
  { id: "r3", status: "running" },
  { id: "r4", status: "ok" },
];

describe("runs selected for a comparison", () => {
  test.each<[string, JobId[], ReturnType<typeof comparePair>]>([
    ["keep the order of the selection, not of the list", ["r4", "r1"], { a: "r4", b: "r1" }],
    ["need a second run", ["r1"], undefined],
    ["need both runs finished", ["r1", "r3"], undefined],
    ["need exactly two runs", ["r1", "r2", "r4"], undefined],
    ["need both runs listed", ["r1", "r9"], undefined],
  ])("%s", (_, selected, expected) => {
    expect(comparePair(selected, RUNS)).toStrictEqual(expected);
  });
});
