import { type FocusSituation, settingFocusStep, type FocusStep } from "../settingFocus";

describe("focus on a setting of All settings", () => {
  test.each<[string, FocusSituation, FocusStep]>([
    ["a shown field takes the focus", { found: true, view: "all", filtered: true }, "focus"],
    ["a closed sheet opens first", { found: false, view: "main", filtered: false }, "open"],
    ["a closed sheet opens also with a search", { found: false, view: "main", filtered: true }, "open"],
    ["a search that hides the field is cleared", { found: false, view: "all", filtered: true }, "reveal"],
    ["an opening sheet leaves the focus to its initial focus", { found: false, view: "all", filtered: false }, "wait"],
  ])("%s", (_name, situation, expected) => {
    expect(settingFocusStep(situation)).toBe(expected);
  });
});
