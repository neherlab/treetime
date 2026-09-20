import { noAlwaysTrueAssertionRule } from "../rules/no-always-true-assertion.ts"
import { ruleTester } from "./rule-tester.ts"

const tester = ruleTester("ts")

tester.run("treetime/no-always-true-assertion", noAlwaysTrueAssertionRule, {
  valid: ["expect(result).toBe(2)", "expect(value).toEqual(expected)"],
  invalid: [
    { code: "expect(1).toBe(1)", errors: [{ messageId: "literal" }] },
    { code: "expect('ok').toBeTruthy()", errors: [{ messageId: "literal" }] },
  ],
})
