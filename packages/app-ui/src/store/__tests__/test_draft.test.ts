import { describe, expect, test } from "vitest";

import { freshDraft, storedDraft } from "../draftSchema";

describe("draft", () => {
  test("a stored draft keeps its fields", () => {
    const draft = { ...freshDraft("clock"), search: "rate", config: { tree: "tree.nwk", clock_rate: 0.003 } };

    expect(storedDraft(draft)).toStrictEqual(draft);
  });

  test("a stored draft whose configuration is not JSON is refused", () => {
    const draft = { ...freshDraft("clock"), config: { tree: undefined } };

    expect(() => storedDraft(draft)).toThrow("Invalid input");
  });
});
