import { noTestResourceAccessRule } from "../rules/no-test-resource-access.ts"
import { ruleTester } from "./rule-tester.ts"

const tester = ruleTester("ts")
const error = { messageId: "resource" }

tester.run("treetime/no-test-resource-access", noTestResourceAccessRule, {
  valid: [
    { code: "const x = 1", filename: "widget.test.ts" },
    { code: "expect(add(1, 2)).toBe(3)", filename: "math.test.ts" },
    { code: "function setTimeout() {} setTimeout()", filename: "widget.test.ts" },
    { code: "const clock = { now() { return 0 } }; clock.now()", filename: "widget.test.ts" },
    { code: "setTimeout(run, 10)", filename: "flow.integration.test.ts" },
    { code: "fetch('/api')", filename: "flow.process.test.ts" },
    { code: "performance.now()", filename: "flow.migration.test.ts" },
    { code: "setInterval(run, 1)", filename: "flow.gm.test.ts" },
    { code: "fetch('/api')", filename: "test_gm_flow.ts" },
  ],
  invalid: [
    { code: "setTimeout(run, 10)", filename: "widget.test.ts", errors: [error] },
    { code: "setInterval(run, 10)", filename: "widget.test.ts", errors: [error] },
    { code: "fetch('/api')", filename: "widget.test.ts", errors: [error] },
    { code: "performance.now()", filename: "widget.test.ts", errors: [error] },
    { code: "Date.now()", filename: "widget.test.ts", errors: [error] },
  ],
})
