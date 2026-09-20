import { defineRule } from "@oxlint/plugins"
import type { ESTree } from "@oxlint/plugins"

import { isFunctionNode } from "./function-shape.ts"
import type { FunctionNode } from "./function-shape.ts"

const LOOP_STATEMENTS = new Set([
  "ForStatement",
  "ForInStatement",
  "ForOfStatement",
  "WhileStatement",
  "DoWhileStatement",
])
const ITERATION_METHODS = new Set([
  "forEach",
  "map",
  "flatMap",
  "filter",
  "some",
  "every",
  "find",
  "reduce",
  "reduceRight",
])

export const noAssertionInLoopRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description:
        "Disallow `expect` and `assert` calls inside a loop or an array iteration callback in a test.",
    },
    messages: {
      inLoop:
        "Assertion inside `{{construct}}`: the run stops at the first failing case and the report does not name it. Use `test.each`, or one test per case.",
    },
  },
  createOnce(context) {
    return {
      CallExpression(node) {
        if (!isAssertion(node)) {
          return
        }
        const construct = enclosingIteration(node)
        if (construct !== undefined) {
          context.report({ node, messageId: "inLoop", data: { construct } })
        }
      },
    }
  },
})

function isAssertion(node: ESTree.CallExpression): boolean {
  const { callee } = node
  if (callee.type === "Identifier") {
    return callee.name === "expect" || callee.name === "assert"
  }
  return (
    callee.type === "MemberExpression" &&
    callee.object.type === "Identifier" &&
    callee.object.name === "assert"
  )
}

function enclosingIteration(node: ESTree.Node): string | undefined {
  let current: ESTree.Node | null = node.parent
  while (current !== null) {
    if (LOOP_STATEMENTS.has(current.type)) {
      return loopKeyword(current.type)
    }
    if (isFunctionNode(current)) {
      const method = iterationMethod(current)
      if (method !== undefined) {
        return `${method}()`
      }
    }
    current = current.parent
  }
  return undefined
}

function iterationMethod(functionNode: FunctionNode): string | undefined {
  const { parent } = functionNode
  if (
    parent.type !== "CallExpression" ||
    !parent.arguments.some((argument) => argument === functionNode) ||
    parent.callee.type !== "MemberExpression" ||
    parent.callee.computed ||
    parent.callee.property.type !== "Identifier" ||
    !ITERATION_METHODS.has(parent.callee.property.name)
  ) {
    return undefined
  }
  return parent.callee.property.name
}

function loopKeyword(type: string): string {
  if (type === "WhileStatement") {
    return "while"
  }
  if (type === "DoWhileStatement") {
    return "do"
  }
  return "for"
}
