import { defineRule } from "@oxlint/plugins";

import { memberRoot } from "./function-shape.ts";

export const preferStrictEqualRule = defineRule({
  meta: {
    type: "suggestion",
    docs: {
      description:
        "Disallow `toEqual` on an expect chain; `toStrictEqual` also checks undefined properties and class identity.",
    },
    messages: {
      toEqual:
        "`toEqual` ignores `undefined` properties, sparse slots, and class identity, so two different values pass. Use `toStrictEqual`.",
    },
  },
  createOnce(context) {
    return {
      CallExpression(node) {
        const { callee } = node;

        if (
          callee.type !== "MemberExpression" ||
          callee.computed ||
          callee.property.type !== "Identifier" ||
          callee.property.name !== "toEqual"
        ) {
          return;
        }

        const root = memberRoot(callee.object);

        if (root.type === "CallExpression" && root.callee.type === "Identifier" && root.callee.name === "expect") {
          context.report({ node: callee.property, messageId: "toEqual" });
        }
      },
    };
  },
});
