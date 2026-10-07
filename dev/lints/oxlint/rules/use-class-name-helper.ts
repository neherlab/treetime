import { defineRule } from "@oxlint/plugins";

export const useClassNameHelperRule = defineRule({
  meta: {
    type: "suggestion",
    docs: {
      description: "Disallow importing `clsx` or `tailwind-merge` directly; use `cn`, the project's class name helper.",
    },
    messages: {
      useCn: "Import `cn`, the project's class name helper, instead of `{{source}}`.",
    },
  },
  createOnce(context) {
    return {
      ImportDeclaration(node) {
        const source = node.source.value;

        if (source === "clsx" || source === "tailwind-merge") {
          context.report({ node, messageId: "useCn", data: { source } });
        }
      },
    };
  },
});
