import { requireSuppressionReasonRule } from "../rules/require-suppression-reason.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

const error = { messageId: "missingReason" };

tester.run("treetime/require-suppression-reason", requireSuppressionReasonRule, {
  valid: [
    "// oxlint-disable-next-line treetime/no-vague-identifiers -- entry name is fixed\nconst utils = 1",
    "// oxlint-disable treetime/foo -- justified\nconst x = 1\n// oxlint-enable treetime/foo",
    "// oxlint-enable treetime/foo\nconst x = 1",
    "// a normal comment\nconst x = 1",
    "const x = 1",
  ],
  invalid: [
    {
      code: "// oxlint-disable-next-line treetime/no-vague-identifiers\nconst utils = 1",
      errors: [error],
    },
    { code: "// oxlint-disable-line treetime/foo\nconst x = 1", errors: [error] },
    {
      code: "// oxlint-disable treetime/foo\nconst x = 1\n// oxlint-enable treetime/foo",
      errors: [error],
    },
    { code: "// oxlint-disable-next-line treetime/foo --\nconst x = 1", errors: [error] },
  ],
});
