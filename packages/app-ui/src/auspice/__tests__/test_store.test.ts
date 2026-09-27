import { CLEAN_START } from "auspice/src/actions/types";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";

import { createAuspiceStore } from "../store";

describe("auspice store", () => {
  beforeEach(() => {
    vi.stubGlobal("window", { innerWidth: 1280, innerHeight: 800, document: { body: { clientHeight: 800 } } });
  });

  afterEach(() => {
    vi.unstubAllGlobals();
  });

  test.each([
    { tips: 4000, skipTreeAnimation: false },
    { tips: 4001, skipTreeAnimation: true },
  ])(
    "a tree of $tips tips sets skipTreeAnimation to $skipTreeAnimation (auspice performanceFlags: above 4000 tips)",
    ({ tips, skipTreeAnimation }) => {
      const store = createAuspiceStore();

      store.dispatch(cleanStart(tips));

      expect(store.getState().controls.performanceFlags).toStrictEqual(
        new Map([["skipTreeAnimation", skipTreeAnimation]]),
      );
    },
  );
});

function cleanStart(fullTipCount: number) {
  return {
    type: CLEAN_START,
    metadata: { loaded: true },
    tree: { loaded: true, nodes: [{ fullTipCount }] },
    controls: {},
    entropy: {},
    measurements: {},
  };
}
