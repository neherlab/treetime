import { defineRule } from "@oxlint/plugins"

export const noModuleLevelMutableRule = defineRule({
  meta: {
    type: "suggestion",
    docs: {
      description:
        "Disallow module-level mutable bindings, which are shared by every importer and every test in the process.",
    },
    messages: {
      moduleMutable:
        "Module-level mutable state is banned. Use a const binding or encapsulate state in a function or store.",
    },
  },
  createOnce(context) {
    return {
      Program(node) {
        for (const statement of node.body) {
          const declaration =
            statement.type === "ExportNamedDeclaration" ? statement.declaration : statement
          if (
            declaration != null &&
            declaration.type === "VariableDeclaration" &&
            declaration.kind !== "const"
          ) {
            context.report({ node: declaration, messageId: "moduleMutable" })
          }
        }
      },
    }
  },
})
