import { noModuleMocksRule } from "../rules/no-module-mocks.ts"
import { ruleTester } from "./rule-tester.ts"

const tester = ruleTester("ts")

tester.run("treetime/no-module-mocks", noModuleMocksRule, {
  valid: ["vi.fn()", "inject(dependency)"],
  invalid: [
    { code: "vi.mock('./service')", errors: [{ messageId: "moduleMock" }] },
    { code: "vi.doMock('./service')", errors: [{ messageId: "moduleMock" }] },
    { code: "vi.unmock('./service')", errors: [{ messageId: "moduleMock" }] },
  ],
})
