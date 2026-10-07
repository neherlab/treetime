import { requireSuppressionReasonRule } from "../rules/require-suppression-reason.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

const error = { messageId: "missingReason" };

tester.run("custom/require-suppression-reason", requireSuppressionReasonRule, {
  valid: [
    "// oxlint-disable-next-line custom/no-vague-identifiers -- entry name is fixed\nconst utils = 1",
    "// oxlint-disable custom/foo -- justified\nconst x = 1\n// oxlint-enable custom/foo",
    "// oxlint-enable custom/foo\nconst x = 1",
    "// a normal comment\nconst x = 1",
    "const x = 1",
  ],
  invalid: [
    {
      code: "// oxlint-disable-next-line custom/no-vague-identifiers\nconst utils = 1",
      errors: [error],
    },
    { code: "// oxlint-disable-line custom/foo\nconst x = 1", errors: [error] },
    {
      code: "// oxlint-disable custom/foo\nconst x = 1\n// oxlint-enable custom/foo",
      errors: [error],
    },
    { code: "// oxlint-disable-next-line custom/foo --\nconst x = 1", errors: [error] },
  ],
});
