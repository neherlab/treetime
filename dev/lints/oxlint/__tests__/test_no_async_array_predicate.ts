import { noAsyncArrayPredicateRule } from "../rules/no-async-array-predicate.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

tester.run("treetime/no-async-array-predicate", noAsyncArrayPredicateRule, {
  valid: ["items.filter((item) => item.ok)", "items.map(async (item) => await load(item))"],
  invalid: [
    {
      code: "items.filter(async (item) => await ok(item))",
      errors: [{ messageId: "asyncPredicate", data: { method: "filter" } }],
    },
    {
      code: "items.some(async (item) => await ok(item))",
      errors: [{ messageId: "asyncPredicate", data: { method: "some" } }],
    },
  ],
});
