import { noLeakyMocksRule } from "../rules/no-leaky-mocks.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

tester.run("custom/no-leaky-mocks", noLeakyMocksRule, {
  valid: ["vi.mock('./service')", "makeStub()"],
  invalid: [
    { code: "vi.spyOn(console, 'log')", errors: [{ messageId: "leakyMock" }] },
    { code: "vi.fn()", errors: [{ messageId: "leakyMock" }] },
    { code: "run(vi.fn<() => void>())", errors: [{ messageId: "leakyMock" }] },
    { code: "expect(() => run(vi.fn())).not.toThrow()", errors: [{ messageId: "leakyMock" }] },
  ],
});
