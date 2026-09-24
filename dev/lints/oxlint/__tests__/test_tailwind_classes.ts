import { tailwindClassesRule } from "../rules/tailwind-classes.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("tsx");

tester.run("web/tailwind-classes", tailwindClassesRule, {
  valid: [
    'const node = <div className="flex" />',
    'const node = <div className="bg-white hover:bg-black" />',
    'const node = <input className="text-black placeholder:text-white" />',
    'const node = <progress className="bg-white [&::-moz-progress-bar]:bg-black" />',
    'const node = <div className="border-white data-[invalid]:border-black" />',
  ],
  invalid: [
    {
      code: 'const node = <div className="totally-not-a-tailwind-class-xyz" />',
      errors: [{ messageId: "unknown", data: { token: "totally-not-a-tailwind-class-xyz" } }],
    },
    {
      code: 'const node = <div className="flex flex" />',
      errors: [{ messageId: "duplicate", data: { token: "flex" } }],
    },
    {
      code: 'const node = <div className="bg-white bg-black" />',
      errors: [{ messageId: "conflict", data: { prior: "bg-white", token: "bg-black" } }],
    },
    {
      code: 'const node = <div className="hover:bg-white hover:bg-black" />',
      errors: [{ messageId: "conflict", data: { prior: "hover:bg-white", token: "hover:bg-black" } }],
    },
  ],
});
