import { noFocusedTestsRule } from "../rules/no-focused-tests.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

tester.run("treetime/no-focused-tests", noFocusedTestsRule, {
  valid: ["test('runs', () => {})", "describe('suite', () => {})"],
  invalid: [
    { code: "fit('focus', () => {})", errors: [{ messageId: "focused" }] },
    { code: "test.only('focus', () => {})", errors: [{ messageId: "focused" }] },
    { code: "describe.only('focus', () => {})", errors: [{ messageId: "focused" }] },
  ],
});
