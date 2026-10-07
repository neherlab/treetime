import { noContractAliasRule } from "../rules/no-contract-alias.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

const OPTIONS = [{ package: "@neherlab/app-contracts" }];

tester.run("custom/no-contract-alias", noContractAliasRule, {
  valid: [
    {
      code: 'import type { RunRecord } from "@neherlab/app-contracts";\nexport function title(record: RunRecord) { return record.title; }',
      options: OPTIONS,
    },
    {
      code: 'import type { RunRecord } from "@neherlab/app-contracts";\ntype ConfigOf = RunRecord["config"];',
      options: OPTIONS,
    },
    {
      code: 'import type { RunRecord } from "@neherlab/app-contracts";\ntype Picked = Pick<RunRecord, "id">;',
      options: OPTIONS,
    },
    { code: 'import type { Local } from "./local";\ntype Other = Local;', options: OPTIONS },
    { code: 'export type { RunRecord } from "@neherlab/app-contracts";', options: OPTIONS },
  ],
  invalid: [
    {
      code: 'import type { RunFile } from "@neherlab/app-contracts";\nexport type RunFileEntry = RunFile;',
      options: OPTIONS,
      errors: [{ messageId: "alias", data: { alias: "RunFileEntry", name: "RunFile" } }],
    },
    {
      code: 'import type { SettingSpec } from "@neherlab/app-contracts";\ntype Spec = Omit<SettingSpec, "examples"> & { examples: string[] };',
      options: OPTIONS,
      errors: [{ messageId: "patched", data: { name: "SettingSpec" } }],
    },
    {
      code: 'export type { ClockResults as ClockData } from "@neherlab/app-contracts";',
      options: OPTIONS,
      errors: [{ messageId: "renamedExport", data: { name: "ClockResults", alias: "ClockData" } }],
    },
    {
      code: 'import type { ClockResults } from "@neherlab/app-contracts";\nexport type { ClockResults as ClockData };',
      options: OPTIONS,
      errors: [{ messageId: "renamedExport", data: { name: "ClockResults", alias: "ClockData" } }],
    },
  ],
});
