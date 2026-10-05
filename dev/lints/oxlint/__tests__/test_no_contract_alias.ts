import { noContractAliasRule } from "../rules/no-contract-alias.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

tester.run("treetime/no-contract-alias", noContractAliasRule, {
  valid: [
    'import type { RunRecord } from "@neherlab/app-contracts";\nexport function title(record: RunRecord) { return record.title; }',
    'import type { RunRecord } from "@neherlab/app-contracts";\ntype ConfigOf = RunRecord["config"];',
    'import type { RunRecord } from "@neherlab/app-contracts";\ntype Picked = Pick<RunRecord, "id">;',
    'import type { Local } from "./local";\ntype Other = Local;',
    'export type { RunRecord } from "@neherlab/app-contracts";',
  ],
  invalid: [
    {
      code: 'import type { RunFile } from "@neherlab/app-contracts";\nexport type RunFileEntry = RunFile;',
      errors: [{ messageId: "alias", data: { alias: "RunFileEntry", name: "RunFile" } }],
    },
    {
      code: 'import type { SettingSpec } from "@neherlab/app-contracts";\ntype Spec = Omit<SettingSpec, "examples"> & { examples: string[] };',
      errors: [{ messageId: "patched", data: { name: "SettingSpec" } }],
    },
    {
      code: 'export type { ClockResults as ClockData } from "@neherlab/app-contracts";',
      errors: [{ messageId: "renamedExport", data: { name: "ClockResults", alias: "ClockData" } }],
    },
    {
      code: 'import type { ClockResults } from "@neherlab/app-contracts";\nexport type { ClockResults as ClockData };',
      errors: [{ messageId: "renamedExport", data: { name: "ClockResults", alias: "ClockData" } }],
    },
  ],
});
