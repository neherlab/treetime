import { noUppercaseTestTitleRule } from "../rules/no-uppercase-test-title.ts"
import { ruleTester } from "./rule-tester.ts"

const tester = ruleTester("ts")

tester.run("treetime/no-uppercase-test-title", noUppercaseTestTitleRule, {
  valid: ["test('adds two numbers', () => {})", "describe('parser', () => {})"],
  invalid: [
    { code: "test('Adds two numbers', () => {})", errors: [{ messageId: "uppercase" }] },
    { code: "describe('Parser', () => {})", errors: [{ messageId: "uppercase" }] },
  ],
})
