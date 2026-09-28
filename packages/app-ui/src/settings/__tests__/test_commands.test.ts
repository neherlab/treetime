import { describe, expect, test } from "vitest";

import { commandSwitchNote } from "../commands";

describe("command switch note", () => {
  test.each([
    { name: "same command", requested: "timetree", loaded: "timetree", expected: undefined },
    {
      name: "other command",
      requested: "timetree",
      loaded: "clock",
      expected: "Switched to Clock signal",
    },
  ] as const)("$name", ({ requested, loaded, expected }) => {
    expect(commandSwitchNote(requested, loaded)).toStrictEqual(expected);
  });
});
