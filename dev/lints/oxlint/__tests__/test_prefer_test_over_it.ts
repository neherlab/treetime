import { preferTestOverItRule } from "../rules/prefer-test-over-it.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

tester.run("treetime/prefer-test-over-it", preferTestOverItRule, {
  valid: ["test('runs', () => {})", "describe('suite', () => {})"],
  invalid: [{ code: "it('runs', () => {})", errors: [{ messageId: "useTest" }] }],
});
