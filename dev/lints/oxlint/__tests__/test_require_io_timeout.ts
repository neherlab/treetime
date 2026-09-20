import { requireIoTimeoutRule } from "../rules/require-io-timeout.ts"
import { ruleTester } from "./rule-tester.ts"

const tester = ruleTester("ts")

tester.run("treetime/require-io-timeout", requireIoTimeoutRule, {
  valid: [
    "fetch(url, { signal: controller.signal })",
    "notFetch(url)",
  ],
  invalid: [
    { code: "fetch(url)", errors: [{ messageId: "needSignal" }] },
    { code: "fetch(url, { method: 'GET' })", errors: [{ messageId: "needSignal" }] },
  ],
})
