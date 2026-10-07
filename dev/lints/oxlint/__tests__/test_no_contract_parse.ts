import { noContractParseRule } from "../rules/no-contract-parse.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

const OPTIONS = [{ package: "@neherlab/app-contracts" }];

tester.run("custom/no-contract-parse", noContractParseRule, {
  valid: [
    { code: 'import { zRun } from "./schemas";\nzRun.parse(value);', options: OPTIONS },
    { code: 'import * as z from "zod";\nz.string().parse(value);', options: OPTIONS },
    { code: 'import { zRun } from "@neherlab/app-contracts";\nconst shape = zRun.shape;', options: OPTIONS },
    {
      code: 'import { zRun } from "@neherlab/app-contracts";\nfunction read(zRun: { parse(v: string): string }) { return zRun.parse("a"); }',
      options: OPTIONS,
    },
  ],
  invalid: [
    {
      code: 'import { zRun } from "@neherlab/app-contracts";\nzRun.parse(value);',
      options: OPTIONS,
      errors: [{ messageId: "contractParse", data: { schema: "zRun", method: "parse" } }],
    },
    {
      code: 'import { zRunEvent as events } from "@neherlab/app-contracts";\nevents.safeParse(value);',
      options: OPTIONS,
      errors: [{ messageId: "contractParse", data: { schema: "events", method: "safeParse" } }],
    },
  ],
});
