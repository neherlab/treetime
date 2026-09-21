import { preferStrictEqualRule } from "../rules/prefer-strict-equal.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

tester.run("treetime/prefer-strict-equal", preferStrictEqualRule, {
  valid: ["expect(value).toStrictEqual(expected)", "queue.toEqual(other)"],
  invalid: [
    { code: "expect(value).toEqual(expected)", errors: [{ messageId: "toEqual" }] },
    { code: "expect(value).resolves.toEqual(expected)", errors: [{ messageId: "toEqual" }] },
  ],
});
