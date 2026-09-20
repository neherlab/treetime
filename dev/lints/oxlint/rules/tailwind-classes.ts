import { readFileSync } from "node:fs"
import { createRequire } from "node:module"
import { dirname, resolve } from "node:path"
import { fileURLToPath, pathToFileURL } from "node:url"

import { defineRule } from "@oxlint/plugins"

interface DesignSystem {
  candidatesToCss(classes: string[]): (string | null)[]
}

interface Stylesheet {
  base: string
  content: string
}

interface TailwindModule {
  __unstable__loadDesignSystem(
    css: string,
    options: { base: string; loadStylesheet(id: string, base: string): Stylesheet },
  ): DesignSystem
}

const PROJECT_ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..", "..", "..")
const CSS_ENTRY = resolve(PROJECT_ROOT, "packages/app-web/src/index.css")

const designSystem = await loadDesignSystem()
const classKeyCache = new Map<string, string | null>()

export const tailwindClassesRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description: "Reject unknown, duplicate, and conflicting Tailwind classes against the real stylesheet.",
    },
    messages: {
      duplicate: "Duplicate Tailwind class `{{token}}`.",
      unknown: "Unknown Tailwind class `{{token}}` is not in the stylesheet.",
      conflict: "Conflicting Tailwind classes `{{prior}}` and `{{token}}` set the same property.",
    },
  },
  createOnce(context) {
    function check(node: { type: string }, value: string): void {
      const tokens = value.split(/\s+/).filter(Boolean)
      const seen = new Set<string>()
      const keyToClass = new Map<string, string>()
      for (const token of tokens) {
        if (seen.has(token)) {
          context.report({ node, messageId: "duplicate", data: { token } })
          continue
        }
        seen.add(token)
        const key = classRuleKey(token)
        if (key === null) {
          context.report({ node, messageId: "unknown", data: { token } })
          continue
        }
        const prior = keyToClass.get(key)
        if (prior !== undefined) {
          context.report({ node, messageId: "conflict", data: { prior, token } })
        } else {
          keyToClass.set(key, token)
        }
      }
    }
    function checkNode(node: { type: string; value?: unknown }): void {
      if (node.type === "Literal" && typeof node.value === "string") {
        check(node, node.value)
      }
    }
    return {
      JSXAttribute(node) {
        if (node.name.type !== "JSXIdentifier" || node.name.name !== "className" || node.value == null) {
          return
        }
        if (node.value.type === "Literal") {
          checkNode(node.value)
        } else if (node.value.type === "JSXExpressionContainer") {
          checkNode(node.value.expression)
        }
      },
      CallExpression(node) {
        if (node.callee.type === "Identifier" && node.callee.name === "cn") {
          for (const argument of node.arguments) {
            checkNode(argument)
          }
        }
      },
    }
  },
})

async function loadDesignSystem(): Promise<DesignSystem> {
  const req = createRequire(CSS_ENTRY)
  const imported: unknown = await import(pathToFileURL(req.resolve("tailwindcss")).href)
  const tailwind = pickTailwindModule(imported)
  const css = readFileSync(CSS_ENTRY, "utf8")
  const loadStylesheet = (id: string, base: string): Stylesheet => {
    let file: string
    if (id === "tailwindcss") {
      file = req.resolve("tailwindcss/index.css")
    } else if (id.startsWith("tailwindcss/")) {
      file = req.resolve(id)
    } else {
      file = resolve(base, id)
    }
    return { base: dirname(file), content: readFileSync(file, "utf8") }
  }
  return tailwind.__unstable__loadDesignSystem(css, { base: dirname(CSS_ENTRY), loadStylesheet })
}

function pickTailwindModule(imported: unknown): TailwindModule {
  const candidate = hasDefault(imported) ? imported.default : imported
  if (!hasLoader(candidate)) {
    throw new Error("tailwindcss does not expose __unstable__loadDesignSystem")
  }
  return candidate
}

function hasDefault(value: unknown): value is { default: unknown } {
  return typeof value === "object" && value !== null && "default" in value
}

function hasLoader(value: unknown): value is TailwindModule {
  return (
    typeof value === "object" &&
    value !== null &&
    typeof Reflect.get(value, "__unstable__loadDesignSystem") === "function"
  )
}

function classRuleKey(className: string): string | null {
  const cached = classKeyCache.get(className)
  if (cached !== undefined) {
    return cached
  }
  const [css] = designSystem.candidatesToCss([className])
  if (css == null) {
    classKeyCache.set(className, null)
    return null
  }
  const context: string[] = []
  const byContext = new Map<string, string[]>()
  for (const raw of css.split("\n")) {
    const line = raw.trim()
    if (line.length === 0) {
      continue
    }
    if (line.endsWith("{")) {
      context.push(line.slice(0, -1).trim())
    } else if (line === "}") {
      context.pop()
    } else {
      const match = /^(--[\w-]+|[a-zA-Z-]+)\s*:/.exec(line)
      if (match) {
        const path = context.slice(1).join(" >> ")
        const props = byContext.get(path) ?? []
        props.push(match[1])
        byContext.set(path, props)
      }
    }
  }
  const key = JSON.stringify(
    [...byContext.entries()].sort().map(([path, props]) => [path, [...props].sort()]),
  )
  classKeyCache.set(className, key)
  return key
}
