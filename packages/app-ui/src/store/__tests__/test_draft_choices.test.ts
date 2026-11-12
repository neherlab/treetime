import { withPickedChoice } from "../draftSchema";

describe("picked choices of a draft", () => {
  test("a pick replaces the earlier pick of the same choice and keeps the other choices", () => {
    const picked = withPickedChoice(
      [
        { choice: "clock-rate", option: "fixed" },
        { choice: "root", option: "reroot" },
      ],
      { choice: "root", option: "outgroup" },
    );

    expect(picked).toStrictEqual([
      { choice: "clock-rate", option: "fixed" },
      { choice: "root", option: "outgroup" },
    ]);
  });

  test("a draft saved without picks takes the first pick", () => {
    expect(withPickedChoice(undefined, { choice: "root", option: "keep" })).toStrictEqual([
      { choice: "root", option: "keep" },
    ]);
  });
});
