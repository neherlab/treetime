import { defineRule } from "@oxlint/plugins";

import { bindingVariable } from "./binding.ts";

const EXEMPT_TEST_FILE = /(?:\.integration|\.process|\.migration|\.gm)\.test\.tsx?$|(?:^|[/\\])test_gm_/u;

const RESOURCE_GLOBALS = new Set([
  "setTimeout",
  "setInterval",
  "setImmediate",
  "queueMicrotask",
  "requestAnimationFrame",
  "fetch",
]);

const CLOCK_OWNERS = new Set(["performance", "Date"]);

export const noTestResourceAccessRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description:
        "Disallow timer, wall-clock, and network access in unit tests. Integration, process, migration, and golden-master tests are exempt.",
    },
    messages: {
      resource:
        "A unit test must not reach a timer, the wall clock, or the network. Inject the value, or move this into an integration, process, migration, or golden-master test.",
    },
  },
  createOnce(context) {
    return {
      CallExpression(node) {
        if (EXEMPT_TEST_FILE.test(context.filename)) {
          return;
        }

        const { callee } = node;

        if (callee.type === "Identifier") {
          if (RESOURCE_GLOBALS.has(callee.name) && bindingVariable(callee, context.sourceCode) === undefined) {
            context.report({ node, messageId: "resource" });
          }

          return;
        }

        if (
          callee.type === "MemberExpression" &&
          !callee.computed &&
          callee.property.type === "Identifier" &&
          callee.property.name === "now" &&
          callee.object.type === "Identifier" &&
          CLOCK_OWNERS.has(callee.object.name) &&
          bindingVariable(callee.object, context.sourceCode) === undefined
        ) {
          context.report({ node, messageId: "resource" });
        }
      },
    };
  },
});
