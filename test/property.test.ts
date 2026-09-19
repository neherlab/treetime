import { fc, test } from "@fast-check/vitest";
import { expect } from "vitest";

test("reversing an array twice yields the original", () => {
  fc.assert(
    fc.property(fc.array(fc.integer()), (xs) => {
      expect(xs.toReversed().toReversed()).toStrictEqual(xs);
    }),
  );
});

test.prop([fc.string()])("splitting on the empty string then joining round-trips", (s) => {
  expect(s.split("").join("")).toBe(s);
});
