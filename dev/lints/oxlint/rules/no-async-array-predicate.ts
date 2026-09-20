import { defineRule } from "@oxlint/plugins"

const PREDICATE_METHODS = new Set([
  "filter",
  "some",
  "every",
  "find",
  "findIndex",
  "findLast",
  "findLastIndex",
  "sort",
  "toSorted",
])

export const noAsyncArrayPredicateRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description:
        "Disallow an async callback in array predicate and comparator methods, where its promise is read as a boolean or a number.",
    },
    messages: {
      asyncPredicate:
        "`{{method}}` reads the callback result as a value, but an async callback returns a promise, and a promise is always truthy. Resolve the values first, then call `{{method}}` with a synchronous callback.",
    },
  },
  createOnce(context) {
    return {
      CallExpression(node) {
        const { callee } = node
        if (
          callee.type !== "MemberExpression" ||
          callee.computed ||
          callee.property.type !== "Identifier" ||
          !PREDICATE_METHODS.has(callee.property.name)
        ) {
          return
        }
        const callback = node.arguments[0]
        if (
          callback != null &&
          (callback.type === "ArrowFunctionExpression" || callback.type === "FunctionExpression") &&
          callback.async
        ) {
          context.report({
            node: callback,
            messageId: "asyncPredicate",
            data: { method: callee.property.name },
          })
        }
      },
    }
  },
})
