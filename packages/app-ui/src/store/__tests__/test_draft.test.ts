import { describe, expect, test } from "vitest";

import { freshDraft, restoredDraft } from "../draftSchema";

describe("draft", () => {
  test("a stored draft replaces the fields it holds", () => {
    const current = freshDraft("timetree");

    expect(
      restoredDraft({ command: "clock", search: "stored", title: "Run name of an older draft" }, current),
    ).toStrictEqual({
      ...current,
      command: "clock",
      search: "stored",
    });
  });

  test("a stored draft with an invalid field is ignored", () => {
    const current = freshDraft("timetree");

    expect(restoredDraft({ search: "stored", view: "sideways" }, current)).toStrictEqual(current);
  });
});
