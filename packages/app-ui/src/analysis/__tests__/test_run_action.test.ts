import type { RunCheck } from "@neherlab/app-contracts";

import { runBlock, SingleRun } from "../runAction";

const ADVICE: RunCheck = { id: "a", level: "advice", text: "A hint", settings: [] };

const BLOCK: RunCheck = { id: "b", level: "block", text: "Sequence length is required", settings: [] };

describe("what keeps the run from starting", () => {
  test.each<[string, Parameters<typeof runBlock>[0], ReturnType<typeof runBlock>]>([
    ["ready", { starting: false, checks: [ADVICE], checkPending: false, valid: true }, undefined],
    ["a start in flight", { starting: true, checks: [ADVICE], checkPending: false, valid: true }, "starting"],
    ["no check yet", { starting: false, checks: undefined, checkPending: false, valid: true }, "checking"],
    ["a check of an earlier edit", { starting: false, checks: [], checkPending: true, valid: true }, "checking"],
    ["a blocking check", { starting: false, checks: [BLOCK], checkPending: false, valid: true }, "blocked"],
    ["an invalid form", { starting: false, checks: [], checkPending: false, valid: false }, "invalid"],
  ])("%s", (_name, readiness, expected) => {
    expect(runBlock(readiness)).toBe(expected);
  });
});

describe("one start at a time", () => {
  test("a second run while the first starts does nothing, and a run after it starts again", async () => {
    const single = new SingleRun();
    const { promise, resolve } = Promise.withResolvers<undefined>();
    let starts = 0;

    const start = async () => {
      starts += 1;
      await promise;
    };

    const first = single.run(start);
    const second = await single.run(start);
    resolve(undefined);
    const firstStarted = await first;

    const third = await single.run(() => {
      starts += 1;

      return Promise.resolve();
    });

    expect([firstStarted, second, third, starts]).toStrictEqual([true, false, true, 2]);
  });
});
