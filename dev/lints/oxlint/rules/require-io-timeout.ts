import { defineRule } from "@oxlint/plugins";

import { calleeName } from "./ast.ts";

export const requireIoTimeoutRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description: "Require `fetch` to pass an AbortSignal through the request options.",
    },
    messages: {
      needSignal: "Network IO must pass an AbortSignal timeout (options.signal).",
    },
  },
  createOnce(context) {
    return {
      CallExpression(node) {
        if (calleeName(node) !== "fetch") {
          return;
        }

        const options = node.arguments[1];

        const hasSignal =
          options != null &&
          options.type === "ObjectExpression" &&
          options.properties.some(
            (property) =>
              property.type === "Property" && property.key.type === "Identifier" && property.key.name === "signal",
          );

        if (!hasSignal) {
          context.report({ node, messageId: "needSignal" });
        }
      },
    };
  },
});
