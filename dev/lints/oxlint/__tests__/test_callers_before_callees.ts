import { callersBeforeCalleesRule } from "../rules/callers-before-callees.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

const error = { messageId: "outOfOrder" };

tester.run("treetime/callers-before-callees", callersBeforeCalleesRule, {
  valid: [
    "function a() { b() }\nfunction b() {}",
    "function a() { const f = () => b(); f() }\nfunction b() {}",
    "function only() { return 1 }",
    "function b() {}\nexport function a() { b() }",
    "function a() { b() }\nfunction b() { a() }",
  ],
  invalid: [
    { code: "function b() {}\nfunction a() { b() }", errors: [error] },
    {
      code: "function c() {}\nfunction b() { c() }\nfunction a() { b() }",
      errors: [error, error],
    },
  ],
});
