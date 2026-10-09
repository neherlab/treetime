import { noFakeSuccessRule } from "../rules/no-fake-success.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

tester.run("custom/no-fake-success", noFakeSuccessRule, {
  valid: [
    "test('adds', () => { expect(add(1, 1)).toBe(2) })",
    "helper(() => {})",
    "test('left to vitest/expect-expect', () => {})",
  ],
  invalid: [
    { code: "describe('nothing', () => {})", errors: [{ messageId: "empty" }] },
    { code: "suite.each([1])('nothing %s', function () {})", errors: [{ messageId: "empty" }] },
  ],
});
