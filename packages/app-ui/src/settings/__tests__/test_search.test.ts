import { commandSettings } from "../catalog";
import { defaultConfig } from "../config";
import { setAt } from "../json";
import { matchingSpecs } from "../search";

const specs = commandSettings("clock").settings;

describe("setting search", () => {
  test("matches the flag of a setting", () => {
    expect(
      matchingSpecs(specs, defaultConfig(specs), "--allow-negative", false).map((spec) => spec.key),
    ).toStrictEqual(["allow_negative_rate"]);
  });

  test("changed only keeps the changed settings", () => {
    const config = setAt(defaultConfig(specs), ["keep_root"], true);

    expect(matchingSpecs(specs, config, "", true).map((spec) => spec.key)).toStrictEqual(["keep_root"]);
  });

  test("every word must match", () => {
    expect(matchingSpecs(specs, defaultConfig(specs), "keep root nonsense", false)).toStrictEqual([]);
  });
});
