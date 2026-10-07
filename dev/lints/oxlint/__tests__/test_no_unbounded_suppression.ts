import { noUnboundedSuppressionRule } from "../rules/no-unbounded-suppression.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

const error = { messageId: "unbounded" };

tester.run("custom/no-unbounded-suppression", noUnboundedSuppressionRule, {
  valid: [
    "// oxlint-disable custom/foo -- reason\nconst x = 1\n// oxlint-enable custom/foo",
    "// oxlint-disable-next-line custom/foo -- reason\nconst x = 1",
    "// oxlint-disable-line custom/foo -- reason\nconst x = 1",
    "// a normal comment\nconst x = 1",
    "const x = 1",
  ],
  invalid: [
    { code: "// oxlint-disable custom/foo -- reason\nconst x = 1", errors: [error] },
    {
      code: "// oxlint-enable custom/foo\n// oxlint-disable custom/foo -- reason\nconst x = 1",
      errors: [error],
    },
  ],
});
