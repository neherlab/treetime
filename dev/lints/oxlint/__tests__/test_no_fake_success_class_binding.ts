import { noFakeSuccessRule } from "../rules/no-fake-success.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts", { Promise: "readonly" });
const swallowed = { messageId: "swallowedError" };

tester.run("no_fake_success_class_binding", noFakeSuccessRule, {
  valid: [
    `Promise.resolve().catch(error => {
  const Holder = class {
    [log(error)]() {}
  }
  return []
})`,
    `Promise.resolve().catch(error => {
  try {
    consume(pending, log(error))
    var pending
  } finally {
    return []
  }
})`,
    `function run() {
  try {
    (0, null)?.method(load())
  } catch {
    return []
  }
}`,
  ],
  invalid: [
    {
      code: `Promise.resolve().catch(error => {
  try {
    class Holder extends load() {
      [log(error)]() {}
    }
  } finally {
    return []
  }
})`,
      errors: [{ ...swallowed, line: 7, column: 4, endLine: 7, endColumn: 13 }],
    },
    {
      code: `Promise.resolve().catch(error => {
  class Holder {
    handle() {
      log(error)
    }
  }
  return []
})`,
      errors: [{ ...swallowed, line: 7, column: 2, endLine: 7, endColumn: 11 }],
    },
    {
      code: `Promise.resolve().catch(error => {
  try {
    switch (2) {
      case 1:
        let handle
      case 2:
        handle(log(error))
    }
  } finally {
    return []
  }
})`,
      errors: [{ ...swallowed, line: 10, column: 4, endLine: 10, endColumn: 13 }],
    },
  ],
});
