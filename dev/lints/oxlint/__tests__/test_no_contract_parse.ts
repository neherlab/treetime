import { noContractParseRule } from "../rules/no-contract-parse.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

tester.run("treetime/no-contract-parse", noContractParseRule, {
  valid: [
    'import { zRun } from "./schemas";\nzRun.parse(value);',
    'import * as z from "zod";\nz.string().parse(value);',
    'import { zRun } from "@neherlab/app-contracts";\nconst shape = zRun.shape;',
    'import { zRun } from "@neherlab/app-contracts";\nfunction read(zRun: { parse(v: string): string }) { return zRun.parse("a"); }',
  ],
  invalid: [
    {
      code: 'import { zRun } from "@neherlab/app-contracts";\nzRun.parse(value);',
      errors: [{ messageId: "contractParse", data: { schema: "zRun", method: "parse" } }],
    },
    {
      code: 'import { zRunEvent as events } from "@neherlab/app-contracts";\nevents.safeParse(value);',
      errors: [{ messageId: "contractParse", data: { schema: "events", method: "safeParse" } }],
    },
  ],
});
