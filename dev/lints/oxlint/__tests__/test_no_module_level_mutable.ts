import { noModuleLevelMutableRule } from "../rules/no-module-level-mutable.ts"
import { ruleTester } from "./rule-tester.ts"

const tester = ruleTester("ts")

tester.run("treetime/no-module-level-mutable", noModuleLevelMutableRule, {
  valid: [
    "const count = 0",
    "export const config = { retries: 3 }",
    "function next() { let count = 0; return count }",
  ],
  invalid: [
    { code: "let count = 0", errors: [{ messageId: "moduleMutable" }] },
    { code: "export let current = undefined", errors: [{ messageId: "moduleMutable" }] },
    { code: "var total = 1", errors: [{ messageId: "moduleMutable" }] },
  ],
})
