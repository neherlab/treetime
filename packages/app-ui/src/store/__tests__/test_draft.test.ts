import { describe, expect, test } from "vitest";

import { freshDraft, restoredDraft } from "../draftSchema";

describe("draft", () => {
  test("a stored draft replaces the fields it holds", () => {
    const current = freshDraft("timetree");

    expect(restoredDraft({ command: "clock", title: "stored", unknown: 1 }, current)).toStrictEqual({
      ...current,
      command: "clock",
      title: "stored",
    });
  });

  test("a stored draft with an invalid field is ignored", () => {
    const current = freshDraft("timetree");

    expect(restoredDraft({ title: "stored", view: "sideways" }, current)).toStrictEqual(current);
  });
});
