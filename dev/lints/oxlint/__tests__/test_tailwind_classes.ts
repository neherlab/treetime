import { tailwindClassesRule } from "../rules/tailwind-classes.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("tsx");

tester.run("web/tailwind-classes", tailwindClassesRule, {
  valid: ['const node = <div className="flex" />'],
  invalid: [
    {
      code: 'const node = <div className="totally-not-a-tailwind-class-xyz" />',
      errors: [{ messageId: "unknown", data: { token: "totally-not-a-tailwind-class-xyz" } }],
    },
    {
      code: 'const node = <div className="flex flex" />',
      errors: [{ messageId: "duplicate", data: { token: "flex" } }],
    },
  ],
});
