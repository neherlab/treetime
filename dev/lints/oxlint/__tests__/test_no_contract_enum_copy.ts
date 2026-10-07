import { noContractEnumCopyRule } from "../rules/no-contract-enum-copy.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

const ENUMS = [
  {
    enums: [
      ["light", "dark", "system"],
      ["valid", "invalid"],
    ],
  },
];

tester.run("custom/no-contract-enum-copy", noContractEnumCopyRule, {
  valid: [
    'status === "valid"',
    { code: 'type Mode = "light" | "custom";', options: ENUMS },
    { code: 'const modes = ["light", "custom"] as const;', options: ENUMS },
    { code: 'const modes = ["light", "dark"];', options: ENUMS },
    { code: 'type One = "light";', options: ENUMS },
    { code: 'a === "light" || b === "dark"', options: ENUMS },
    { code: 'theme === "light"', options: ENUMS },
  ],
  invalid: [
    {
      code: 'type Theme = "light" | "dark";',
      options: ENUMS,
      errors: [{ messageId: "enumCopy", data: { values: "light, dark" } }],
    },
    {
      code: 'const themes = ["light", "dark", "system"] as const;',
      options: ENUMS,
      errors: [{ messageId: "enumCopy", data: { values: "light, dark, system" } }],
    },
    {
      code: 'const known = status === "valid" || status === "invalid";',
      options: ENUMS,
      errors: [{ messageId: "enumCopy", data: { values: "valid, invalid" } }],
    },
    {
      code: 'const other = theme !== "light" && "dark" !== theme;',
      options: ENUMS,
      errors: [{ messageId: "enumCopy", data: { values: "light, dark" } }],
    },
    {
      code: 'type Status = "valid" | "invalid";',
      errors: [{ messageId: "enumCopy", data: { values: "valid, invalid" } }],
    },
  ],
});
