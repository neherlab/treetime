import { defineRule } from "@oxlint/plugins";
import * as z from "zod";

const CN_PACKAGE = "cn";

const zOptions = z.object({ module: z.string() });

export const useThemedCnRule = defineRule({
  meta: {
    type: "suggestion",
    docs: {
      description: "Require importing `cn` from the themed ui module instead of the `cn` package.",
    },
    messages: {
      themedCn: "Import `cn` from the themed ui/cn module, not from the `cn` package directly.",
    },
    schema: [
      {
        type: "object",
        properties: { module: { type: "string" } },
        required: ["module"],
        additionalProperties: false,
      },
    ],
  },
  createOnce(context) {
    return {
      ImportDeclaration(node) {
        const options = zOptions.safeParse(context.options[0]);
        const isThemedModule = options.success && context.filename.replaceAll("\\", "/").endsWith(options.data.module);

        if (isThemedModule) {
          return;
        }

        if (node.source.value === CN_PACKAGE) {
          context.report({ node, messageId: "themedCn" });
        }
      },
    };
  },
});
