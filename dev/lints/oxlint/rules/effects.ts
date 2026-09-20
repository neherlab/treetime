import type { ESTree } from "@oxlint/plugins"

export const SUBSCRIPTION_CALLEES = new Set([
  "addEventListener",
  "removeEventListener",
  "subscribe",
  "observe",
  "setInterval",
  "setTimeout",
  "requestAnimationFrame",
])

export interface EffectRecord {
  node: ESTree.CallExpression
  callback: FunctionNode
  deps: ESTree.Expression | ESTree.SpreadElement | undefined
}

export interface EffectState {
  setterNames: Set<string>
  effects: EffectRecord[]
}

type FunctionNode = ESTree.ArrowFunctionExpression | ESTree.FunctionExpression

export function makeEffectState(): {
  state: EffectState
  visitors: {
    VariableDeclarator(node: ESTree.VariableDeclarator): void
    CallExpression(node: ESTree.CallExpression): void
  }
} {
  const state: EffectState = { setterNames: new Set(), effects: [] }
  return {
    state,
    visitors: {
      VariableDeclarator(node) {
        if (node.id.type === "ArrayPattern" && isUseStateCall(node.init)) {
          const setter = node.id.elements[1]
          if (setter != null && setter.type === "Identifier") {
            state.setterNames.add(setter.name)
          }
        }
      },
      CallExpression(node) {
        if (!isEffectCall(node)) {
          return
        }
        const callback = node.arguments[0]
        const deps = node.arguments[1]
        if (
          callback != null &&
          (callback.type === "ArrowFunctionExpression" || callback.type === "FunctionExpression")
        ) {
          state.effects.push({ node, callback, deps })
        }
      },
    },
  }
}

export interface EffectBodyFacts {
  setterCalls: string[]
  hasCleanup: boolean
  statements: ESTree.Statement[]
}

export function effectBodyFacts(callback: FunctionNode, setterNames: Set<string>): EffectBodyFacts {
  const statements = callback.body.type === "BlockStatement" ? callback.body.body : []
  const calleeNames: string[] = []
  collectCalleeNames(callback.body, calleeNames)
  const setterCalls = calleeNames.filter((name) => setterNames.has(name))
  const hasCleanup = statements.some((statement) => statement.type === "ReturnStatement")
  return { setterCalls, hasCleanup, statements }
}

export function hasEmptyDeps(deps: ESTree.Expression | ESTree.SpreadElement | undefined): boolean {
  return deps != null && deps.type === "ArrayExpression" && deps.elements.length === 0
}

export function hasNonEmptyDeps(deps: ESTree.Expression | ESTree.SpreadElement | undefined): boolean {
  return deps != null && deps.type === "ArrayExpression" && deps.elements.length > 0
}

function isUseStateCall(node: ESTree.Expression | null | undefined): boolean {
  if (node == null || node.type !== "CallExpression") {
    return false
  }
  const { callee } = node
  if (callee.type === "Identifier") {
    return callee.name === "useState"
  }
  return (
    callee.type === "MemberExpression" &&
    callee.property.type === "Identifier" &&
    callee.property.name === "useState"
  )
}

function isEffectCall(node: ESTree.CallExpression): boolean {
  const { callee } = node
  if (callee.type === "Identifier") {
    return callee.name === "useEffect"
  }
  return (
    callee.type === "MemberExpression" &&
    callee.property.type === "Identifier" &&
    callee.property.name === "useEffect"
  )
}

function collectCalleeNames(node: unknown, names: string[]): void {
  if (!isNode(node)) {
    return
  }
  if (node.type === "CallExpression" && node.callee.type === "Identifier") {
    names.push(node.callee.name)
  }
  for (const [key, child] of Object.entries(node)) {
    if (key === "parent") {
      continue
    }
    if (Array.isArray(child)) {
      for (const item of child) {
        collectCalleeNames(item, names)
      }
    } else {
      collectCalleeNames(child, names)
    }
  }
}

function isNode(value: unknown): value is ESTree.Node {
  return typeof value === "object" && value !== null && typeof Reflect.get(value, "type") === "string"
}
