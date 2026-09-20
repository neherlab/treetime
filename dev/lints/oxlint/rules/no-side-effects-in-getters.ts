import { defineRule } from "@oxlint/plugins"
import type { ESTree, SourceCode } from "@oxlint/plugins"

import { bindingVariable } from "./binding.ts"

const MUTATORS = new Set([
  "push",
  "pop",
  "shift",
  "unshift",
  "splice",
  "sort",
  "reverse",
  "fill",
  "copyWithin",
  "set",
  "add",
  "delete",
  "clear",
])
const NESTED_FUNCTIONS = new Set([
  "ArrowFunctionExpression",
  "FunctionExpression",
  "FunctionDeclaration",
])

export const noSideEffectsInGettersRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description: "Disallow mutation of outer state inside a getter.",
    },
    messages: {
      sideEffect:
        "A getter must not mutate state outside itself. Return a value and move the mutation into an explicit command.",
    },
  },
  createOnce(context) {
    function checkGetter(fn: ESTree.Node): void {
      if (!("body" in fn) || fn.body === null || fn.body === undefined) {
        return
      }
      const body = fn.body
      if (typeof body !== "object" || Reflect.get(body, "type") !== "BlockStatement") {
        return
      }
      const range: [number, number] = [fn.start, fn.end]
      for (const node of statementsExcludingNestedFunctions(body)) {
        if (node.type === "AssignmentExpression" && isOuterTarget(node.left, range, context.sourceCode)) {
          context.report({ node, messageId: "sideEffect" })
        } else if (
          node.type === "UpdateExpression" &&
          isOuterTarget(node.argument, range, context.sourceCode)
        ) {
          context.report({ node, messageId: "sideEffect" })
        } else if (node.type === "CallExpression" && isMutatingCall(node, range, context.sourceCode)) {
          context.report({ node, messageId: "sideEffect" })
        }
      }
    }
    return {
      MethodDefinition(node) {
        if (node.kind === "get") {
          checkGetter(node.value)
        }
      },
      Property(node) {
        if (node.kind === "get") {
          checkGetter(node.value)
        }
      },
    }
  },
})

function isMutatingCall(
  node: ESTree.CallExpression,
  range: [number, number],
  sourceCode: SourceCode,
): boolean {
  const { callee } = node
  if (callee.type !== "MemberExpression" || callee.computed || callee.property.type !== "Identifier") {
    return false
  }
  return MUTATORS.has(callee.property.name) && isOuterTarget(callee.object, range, sourceCode)
}

function isOuterTarget(target: ESTree.Node, range: [number, number], sourceCode: SourceCode): boolean {
  let current = target
  while (current.type === "MemberExpression") {
    current = current.object
  }
  if (current.type === "ThisExpression") {
    return true
  }
  if (current.type !== "Identifier") {
    return false
  }
  const variable = bindingVariable(current, sourceCode)
  if (variable === undefined) {
    return true
  }
  return !variable.defs.some((def) => def.node.start >= range[0] && def.node.end <= range[1])
}

function statementsExcludingNestedFunctions(root: ESTree.Node): ESTree.Node[] {
  const found: ESTree.Node[] = []
  const stack: ESTree.Node[] = [root]
  while (stack.length > 0) {
    const node = stack.pop()
    if (node === undefined) {
      continue
    }
    found.push(node)
    for (const [key, value] of Object.entries(node)) {
      if (key === "parent") {
        continue
      }
      if (Array.isArray(value)) {
        for (const item of value) {
          if (isNode(item) && !NESTED_FUNCTIONS.has(item.type)) {
            stack.push(item)
          }
        }
      } else if (isNode(value) && !NESTED_FUNCTIONS.has(value.type)) {
        stack.push(value)
      }
    }
  }
  return found
}

function isNode(value: unknown): value is ESTree.Node {
  return typeof value === "object" && value !== null && typeof Reflect.get(value, "type") === "string"
}
