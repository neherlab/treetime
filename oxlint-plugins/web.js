import { readFileSync } from "node:fs";
import { createRequire } from "node:module";
import { dirname, resolve } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";

const PROJECT_ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const CSS_ENTRY = resolve(PROJECT_ROOT, "packages/app-web/src/index.css");

const SUBSCRIPTION_CALLEES = new Set([
  "addEventListener",
  "removeEventListener",
  "subscribe",
  "observe",
  "setInterval",
  "setTimeout",
  "requestAnimationFrame",
]);

const designSystem = await loadDesignSystem();
const classKeyCache = new Map();

async function loadDesignSystem() {
  const req = createRequire(CSS_ENTRY);
  const tailwind = await import(pathToFileURL(req.resolve("tailwindcss")).href);
  const load = (tailwind.default ?? tailwind).__unstable__loadDesignSystem;
  const css = readFileSync(CSS_ENTRY, "utf8");
  const loadStylesheet = (id, base) => {
    let file;
    if (id === "tailwindcss") {
      file = req.resolve("tailwindcss/index.css");
    } else if (id.startsWith("tailwindcss/")) {
      file = req.resolve(id);
    } else {
      file = resolve(base, id);
    }
    return { base: dirname(file), content: readFileSync(file, "utf8") };
  };
  return load(css, { base: dirname(CSS_ENTRY), loadStylesheet });
}

function classRuleKey(className) {
  if (classKeyCache.has(className)) {
    return classKeyCache.get(className);
  }
  const [css] = designSystem.candidatesToCss([className]);
  if (css == null) {
    classKeyCache.set(className, null);
    return null;
  }
  const context = [];
  const byContext = new Map();
  for (const raw of css.split("\n")) {
    const line = raw.trim();
    if (line.length === 0) {
      continue;
    }
    if (line.endsWith("{")) {
      context.push(line.slice(0, -1).trim());
    } else if (line === "}") {
      context.pop();
    } else {
      const match = /^(--[\w-]+|[a-zA-Z-]+)\s*:/.exec(line);
      if (match) {
        const path = context.slice(1).join(" >> ");
        const props = byContext.get(path) ?? [];
        props.push(match[1]);
        byContext.set(path, props);
      }
    }
  }
  const key = JSON.stringify(
    [...byContext.entries()].sort().map(([path, props]) => [path, [...props].sort()]),
  );
  classKeyCache.set(className, key);
  return key;
}

const tailwindClasses = {
  meta: { name: "tailwind-classes" },
  create(context) {
    const check = (node, value) => {
      const tokens = value.split(/\s+/).filter(Boolean);
      const seen = new Set();
      const keyToClass = new Map();
      for (const token of tokens) {
        if (seen.has(token)) {
          context.report({ node, message: `Duplicate Tailwind class "${token}".` });
          continue;
        }
        seen.add(token);
        const key = classRuleKey(token);
        if (key === null) {
          context.report({
            node,
            message: `Unknown Tailwind class "${token}" is not in the stylesheet.`,
          });
          continue;
        }
        const prior = keyToClass.get(key);
        if (prior !== undefined) {
          context.report({
            node,
            message: `Conflicting Tailwind classes "${prior}" and "${token}" set the same property.`,
          });
        } else {
          keyToClass.set(key, token);
        }
      }
    };
    const checkNode = (node) => {
      if (node.type === "Literal" && typeof node.value === "string") {
        check(node, node.value);
      }
    };
    return {
      JSXAttribute(node) {
        if (node.name.type !== "JSXIdentifier" || node.name.name !== "className" || !node.value) {
          return;
        }
        if (node.value.type === "Literal") {
          checkNode(node.value);
        } else if (node.value.type === "JSXExpressionContainer") {
          checkNode(node.value.expression);
        }
      },
      CallExpression(node) {
        if (node.callee.type === "Identifier" && node.callee.name === "cn") {
          for (const arg of node.arguments) {
            checkNode(arg);
          }
        }
      },
    };
  },
};

const CN_SOURCES = new Set(["clsx", "tailwind-merge", "cn"]);

const useThemedCn = {
  meta: { name: "use-themed-cn" },
  create(context) {
    const filename = context.filename ?? context.getFilename?.() ?? "";
    const isThemedModule = filename.replace(/\\/g, "/").endsWith("app-ui/src/ui/cn.ts");
    return {
      ImportDeclaration(node) {
        if (isThemedModule) {
          return;
        }
        if (typeof node.source.value === "string" && CN_SOURCES.has(node.source.value)) {
          context.report({
            node,
            message: `Import "cn" from the themed ui/cn module, not "${node.source.value}" directly.`,
          });
        }
      },
    };
  },
};

function isUseStateCall(node) {
  if (!node || node.type !== "CallExpression") {
    return false;
  }
  const callee = node.callee;
  if (callee.type === "Identifier") {
    return callee.name === "useState";
  }
  return (
    callee.type === "MemberExpression" &&
    callee.property.type === "Identifier" &&
    callee.property.name === "useState"
  );
}

function collectCalleeNames(node, names) {
  if (!node || typeof node.type !== "string") {
    return;
  }
  if (node.type === "CallExpression" && node.callee.type === "Identifier") {
    names.push(node.callee.name);
  }
  for (const key of Object.keys(node)) {
    if (key === "parent") {
      continue;
    }
    const child = node[key];
    if (Array.isArray(child)) {
      for (const item of child) {
        if (item && typeof item.type === "string") {
          collectCalleeNames(item, names);
        }
      }
    } else if (child && typeof child.type === "string") {
      collectCalleeNames(child, names);
    }
  }
}

function effectRules() {
  const setterNames = new Set();
  const effects = [];

  const record = {
    VariableDeclarator(node) {
      if (node.id.type === "ArrayPattern" && isUseStateCall(node.init)) {
        const setter = node.id.elements[1];
        if (setter && setter.type === "Identifier") {
          setterNames.add(setter.name);
        }
      }
    },
    CallExpression(node) {
      const callee = node.callee;
      const isEffect =
        (callee.type === "Identifier" && callee.name === "useEffect") ||
        (callee.type === "MemberExpression" &&
          callee.property.type === "Identifier" &&
          callee.property.name === "useEffect");
      if (!isEffect) {
        return;
      }
      const cb = node.arguments[0];
      const deps = node.arguments[1];
      if (cb && (cb.type === "ArrowFunctionExpression" || cb.type === "FunctionExpression")) {
        effects.push({ node, cb, deps });
      }
    },
  };
  return { setterNames, effects, record };
}

function effectBodyFacts(cb, setterNames) {
  const body = cb.body.type === "BlockStatement" ? cb.body.body : [{ type: "ReturnStatement" }];
  const calleeNames = [];
  collectCalleeNames(cb.body, calleeNames);
  const setterCalls = calleeNames.filter((name) => setterNames.has(name));
  const hasCleanup = body.some((stmt) => stmt.type === "ReturnStatement");
  const statements = cb.body.type === "BlockStatement" ? cb.body.body : [];
  return { setterCalls, hasCleanup, statements };
}

function makeEffectRule(name, evaluate) {
  return {
    meta: { name },
    create(context) {
      const state = effectRules();
      return {
        ...state.record,
        "Program:exit"() {
          for (const effect of state.effects) {
            evaluate(context, effect, state.setterNames);
          }
        },
      };
    },
  };
}

const noStateInEffectInitializer = makeEffectRule(
  "no-state-in-effect-initializer",
  (context, { node, cb, deps }, setterNames) => {
    const emptyDeps = deps && deps.type === "ArrayExpression" && deps.elements.length === 0;
    if (!emptyDeps) {
      return;
    }
    const { setterCalls } = effectBodyFacts(cb, setterNames);
    if (setterCalls.length > 0) {
      context.report({
        node,
        message:
          "Do not initialize state inside a mount effect. Pass the initial value to useState instead.",
      });
    }
  },
);

const noChainedStateUpdates = makeEffectRule(
  "no-chained-state-updates",
  (context, { node, cb, deps }, setterNames) => {
    const hasDeps = deps && deps.type === "ArrayExpression" && deps.elements.length > 0;
    if (!hasDeps) {
      return;
    }
    const { setterCalls, hasCleanup, statements } = effectBodyFacts(cb, setterNames);
    if (hasCleanup || setterCalls.length === 0 || statements.length === 0) {
      return;
    }
    const onlyStateUpdates = statements.every(
      (stmt) =>
        stmt.type === "ExpressionStatement" &&
        stmt.expression.type === "CallExpression" &&
        stmt.expression.callee.type === "Identifier" &&
        setterNames.has(stmt.expression.callee.name),
    );
    if (onlyStateUpdates) {
      context.report({
        node,
        message:
          "Do not chain state updates in an effect. Derive the value during render or update it in the event handler.",
      });
    }
  },
);

const noEventHandlerEffect = makeEffectRule(
  "no-event-handler-effect",
  (context, { node, cb, deps }, setterNames) => {
    const hasDeps = deps && deps.type === "ArrayExpression" && deps.elements.length > 0;
    if (!hasDeps || cb.body.type !== "BlockStatement" || cb.body.body.length !== 1) {
      return;
    }
    const [stmt] = cb.body.body;
    if (
      stmt.type === "ExpressionStatement" &&
      stmt.expression.type === "CallExpression" &&
      stmt.expression.callee.type === "Identifier" &&
      !setterNames.has(stmt.expression.callee.name) &&
      !SUBSCRIPTION_CALLEES.has(stmt.expression.callee.name)
    ) {
      context.report({
        node,
        message:
          "This effect only forwards to an event handler. Call the handler from the event instead of an effect.",
      });
    }
  },
);

export default {
  meta: { name: "web" },
  rules: {
    "tailwind-classes": tailwindClasses,
    "use-themed-cn": useThemedCn,
    "no-state-in-effect-initializer": noStateInEffectInitializer,
    "no-chained-state-updates": noChainedStateUpdates,
    "no-event-handler-effect": noEventHandlerEffect,
  },
};
