import { noAssertionInLoopRule } from "../rules/no-assertion-in-loop.ts"
import { ruleTester } from "./rule-tester.ts"

const tester = ruleTester("ts")

tester.run("treetime/no-assertion-in-loop", noAssertionInLoopRule, {
  valid: [
    "test('adds', () => { expect(add(1, 1)).toBe(2) })",
    "for (const item of items) { collect(item) }",
  ],
  invalid: [
    {
      code: "for (const item of items) { expect(item).toBe(1) }",
      errors: [{ messageId: "inLoop", data: { construct: "for" } }],
    },
    {
      code: "items.forEach((item) => { expect(item).toBe(1) })",
      errors: [{ messageId: "inLoop", data: { construct: "forEach()" } }],
    },
  ],
})
