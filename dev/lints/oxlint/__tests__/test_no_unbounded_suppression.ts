import { noUnboundedSuppressionRule } from "../rules/no-unbounded-suppression.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

const error = { messageId: "unbounded" };

tester.run("treetime/no-unbounded-suppression", noUnboundedSuppressionRule, {
  valid: [
    "// oxlint-disable treetime/foo -- reason\nconst x = 1\n// oxlint-enable treetime/foo",
    "// oxlint-disable-next-line treetime/foo -- reason\nconst x = 1",
    "// oxlint-disable-line treetime/foo -- reason\nconst x = 1",
    "// a normal comment\nconst x = 1",
    "const x = 1",
  ],
  invalid: [
    { code: "// oxlint-disable treetime/foo -- reason\nconst x = 1", errors: [error] },
    {
      code: "// oxlint-enable treetime/foo\n// oxlint-disable treetime/foo -- reason\nconst x = 1",
      errors: [error],
    },
  ],
});
