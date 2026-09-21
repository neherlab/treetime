import { noDisabledTestsRule } from "../rules/no-disabled-tests.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

tester.run("treetime/no-disabled-tests", noDisabledTestsRule, {
  valid: ["test('runs', () => {})", "describe('suite', () => {})"],
  invalid: [
    { code: "xit('skip', () => {})", errors: [{ messageId: "disabled" }] },
    { code: "test.skip('skip', () => {})", errors: [{ messageId: "disabled" }] },
    { code: "describe.skip('skip', () => {})", errors: [{ messageId: "disabled" }] },
  ],
});
