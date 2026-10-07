import { noDeclareGlobalRule } from "../rules/no-declare-global.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

tester.run("custom/no-declare-global", noDeclareGlobalRule, {
  valid: ["declare module 'x' { export const a: number; }", "namespace App { export const a = 1; }", "export {};"],
  invalid: [
    {
      code: "export {};\ndeclare global { interface Window { host: string } }",
      errors: [{ messageId: "declareGlobal" }],
    },
  ],
});
