import type { ESTree } from "@oxlint/plugins"

export type FunctionNode = ESTree.ArrowFunctionExpression | ESTree.Function

export function isFunctionNode(node: ESTree.Node): node is FunctionNode {
  return (
    node.type === "ArrowFunctionExpression" ||
    node.type === "FunctionDeclaration" ||
    node.type === "FunctionExpression"
  )
}

export function memberRoot(expression: ESTree.Expression): ESTree.Expression {
  let current = expression
  while (current.type === "MemberExpression" || current.type === "ChainExpression") {
    current = current.type === "ChainExpression" ? current.expression : current.object
  }
  return current
}
