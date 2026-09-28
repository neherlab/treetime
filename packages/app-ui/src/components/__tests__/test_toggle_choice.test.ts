import { assert, constantFrom, property, subarray } from "fast-check";
import { describe, expect, test } from "vitest";

import { pressedChoice, toggledChoice } from "../toggleChoice";

const OPTIONS = ["all", "warnings", "stages"] as const;

describe("toggled choice", () => {
  test("pressing another option picks it", () => {
    expect(toggledChoice(OPTIONS, "all", ["all", "stages"])).toBe("stages");
  });

  test("releasing the pressed option keeps the choice", () => {
    expect(toggledChoice(OPTIONS, "warnings", [])).toBeUndefined();
  });

  test("a pressed value outside the options is ignored", () => {
    expect(toggledChoice(OPTIONS, "all", ["all", "unknown"])).toBeUndefined();
  });

  test("the picked value is always an option other than the current one that is pressed", () => {
    assert(
      property(constantFrom(...OPTIONS), subarray([...OPTIONS, "unknown"]), (current, pressed) => {
        const picked = toggledChoice(OPTIONS, current, pressed);
        const candidates = OPTIONS.filter((option) => option !== current && pressed.includes(option));

        expect(picked === undefined ? candidates.length === 0 : candidates.includes(picked)).toBe(true);
      }),
    );
  });
});

describe("pressed choice", () => {
  test("the pressed option is picked", () => {
    expect(pressedChoice(OPTIONS, ["stages"])).toBe("stages");
  });

  test("releasing every option picks nothing", () => {
    expect(pressedChoice(OPTIONS, [])).toBeUndefined();
  });

  test("a pressed value outside the options is ignored", () => {
    expect(pressedChoice(OPTIONS, ["unknown"])).toBeUndefined();
  });

  test("the picked value is always a pressed option", () => {
    assert(
      property(subarray([...OPTIONS, "unknown"]), (pressed) => {
        const picked = pressedChoice(OPTIONS, pressed);
        const candidates = OPTIONS.filter((option) => pressed.includes(option));

        expect(picked === undefined ? candidates.length === 0 : candidates.includes(picked)).toBe(true);
      }),
    );
  });
});
