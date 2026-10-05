import { defineRule } from "@oxlint/plugins";

import { bindingImport } from "./binding.ts";
import { isContractsSource } from "./contracts.ts";

const PARSE_METHODS = new Set(["parse", "safeParse", "parseAsync", "safeParseAsync"]);

export const noContractParseRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description:
        "Disallow parsing with a generated schema outside the trust boundaries; the generated client already validates every response.",
    },
    messages: {
      contractParse:
        "`{{schema}}.{{method}}()` parses data that already has its generated type. Use the typed value; parse only where untrusted data enters the app.",
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
          !PARSE_METHODS.has(callee.property.name) ||
          callee.object.type !== "Identifier"
        ) {
          return;
        }

        const binding = bindingImport(callee.object, context.sourceCode);

        if (binding !== undefined && isContractsSource(binding.source)) {
          context.report({
            node,
            messageId: "contractParse",
            data: { schema: callee.object.name, method: callee.property.name },
          });
        }
      },
    };
  },
});
