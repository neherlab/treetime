import { noFakeSuccessRule } from "../rules/no-fake-success.ts"
import { ruleTester } from "./rule-tester.ts"

const tester = ruleTester("ts")

tester.run("treetime/no-fake-success", noFakeSuccessRule, {
  valid: [
    "test('adds', () => { expect(add(1, 1)).toBe(2) })",
    "helper(() => {})",
  ],
  invalid: [
    { code: "test('does nothing', () => {})", errors: [{ messageId: "empty" }] },
    { code: "it('todo', function () {})", errors: [{ messageId: "empty" }] },
  ],
})
