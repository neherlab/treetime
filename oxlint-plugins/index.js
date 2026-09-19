const VERSIONED_NAME = /(?:_v\d+|_old|_new|_fixed|_deprecated|_copy\d*|_backup|_bak|_temp|_tmp)$/i;
const VAGUE_NAMES = new Set(["utils", "helpers", "misc", "tmp"]);
const TYPOGRAPHIC = /[\u2018\u2019\u201C\u201D\u2013\u2014]|[\u{1F000}-\u{1FAFF}]|[\u{2600}-\u{27BF}]|[\u{2B00}-\u{2BFF}]|\uFE0F/u;
const TEST_CALLERS = new Set(["it", "test", "describe", "bench", "suite"]);

function calleeName(node) {
  if (node.callee.type === "Identifier") {
    return node.callee.name;
  }
  if (node.callee.type === "MemberExpression" && node.callee.property.type === "Identifier") {
    return node.callee.property.name;
  }
  return undefined;
}

function rootCalleeName(node) {
  let current = node.callee;
  while (current.type === "MemberExpression") {
    current = current.object;
  }
  return current.type === "Identifier" ? current.name : undefined;
}

function memberChain(node) {
  const parts = [];
  let current = node.callee;
  while (current.type === "MemberExpression") {
    if (current.property.type === "Identifier") {
      parts.unshift(current.property.name);
    }
    current = current.object;
  }
  if (current.type === "Identifier") {
    parts.unshift(current.name);
  }
  return parts;
}

function declaredName(node) {
  const id = node.id ?? node.key;
  return id && id.type === "Identifier" ? id : undefined;
}

const noModuleLevelMutable = {
  meta: { name: "no-module-level-mutable" },
  create(context) {
    return {
      Program(node) {
        for (const statement of node.body) {
          const decl =
            statement.type === "ExportNamedDeclaration" ? statement.declaration : statement;
          if (decl && decl.type === "VariableDeclaration" && decl.kind !== "const") {
            context.report({
              node: decl,
              message:
                "Module-level mutable state is banned. Use a const binding or encapsulate state in a function or store.",
            });
          }
        }
      },
    };
  },
};

const noVagueIdentifiers = {
  meta: { name: "no-vague-identifiers" },
  create(context) {
    function check(node) {
      const id = declaredName(node);
      if (id && VAGUE_NAMES.has(id.name)) {
        context.report({
          node: id,
          message: `"${id.name}" is a grab-bag name. Name the module or binding for its single responsibility.`,
        });
      }
    }
    return {
      VariableDeclarator: check,
      FunctionDeclaration: check,
      ClassDeclaration: check,
      TSInterfaceDeclaration: check,
      TSTypeAliasDeclaration: check,
      TSModuleDeclaration: check,
      ImportDefaultSpecifier(node) {
        if (VAGUE_NAMES.has(node.local.name)) {
          context.report({ node: node.local, message: `"${node.local.name}" is a grab-bag name.` });
        }
      },
      ImportNamespaceSpecifier(node) {
        if (VAGUE_NAMES.has(node.local.name)) {
          context.report({ node: node.local, message: `"${node.local.name}" is a grab-bag name.` });
        }
      },
    };
  },
};

const noVersionedNames = {
  meta: { name: "no-versioned-names" },
  create(context) {
    function check(node) {
      const id = declaredName(node);
      if (id && VERSIONED_NAME.test(id.name)) {
        context.report({
          node: id,
          message: `"${id.name}" carries a version suffix. Names must describe the final state, not edit history.`,
        });
      }
    }
    return {
      VariableDeclarator: check,
      FunctionDeclaration: check,
      ClassDeclaration: check,
      TSInterfaceDeclaration: check,
      TSTypeAliasDeclaration: check,
    };
  },
};

const noTypographicCharacters = {
  meta: { name: "no-typographic-characters" },
  create(context) {
    return {
      Literal(node) {
        if (typeof node.value === "string" && TYPOGRAPHIC.test(node.value)) {
          context.report({
            node,
            message: "Curly quotes, dashes, and emoji are banned. Use straight ASCII punctuation.",
          });
        }
      },
      TemplateElement(node) {
        const raw = node.value.raw;
        if (typeof raw === "string" && TYPOGRAPHIC.test(raw)) {
          context.report({
            node,
            message: "Curly quotes, dashes, and emoji are banned. Use straight ASCII punctuation.",
          });
        }
      },
    };
  },
};

const noExecShellString = {
  meta: { name: "no-exec-shell-string" },
  create(context) {
    return {
      CallExpression(node) {
        const name = calleeName(node);
        if ((name === "exec" || name === "execSync") && node.arguments.length > 0) {
          const first = node.arguments[0];
          if (first.type === "Literal" || first.type === "TemplateLiteral") {
            context.report({
              node,
              message:
                "exec with a shell string is banned. Use execFile or spawn with an argument array.",
            });
          }
        }
      },
    };
  },
};

const noClassNameHelperBypass = {
  meta: { name: "use-class-name-helper" },
  create(context) {
    return {
      ImportDeclaration(node) {
        const source = node.source.value;
        if (source === "clsx" || source === "tailwind-merge") {
          context.report({
            node,
            message: `Import the themed cn helper instead of "${source}" directly.`,
          });
        }
      },
    };
  },
};

const requireIoTimeout = {
  meta: { name: "require-io-timeout" },
  create(context) {
    return {
      CallExpression(node) {
        if (calleeName(node) !== "fetch") {
          return;
        }
        const options = node.arguments[1];
        const hasSignal =
          options &&
          options.type === "ObjectExpression" &&
          options.properties.some(
            (property) =>
              property.type === "Property" &&
              property.key.type === "Identifier" &&
              property.key.name === "signal",
          );
        if (!hasSignal) {
          context.report({
            node,
            message: "Network IO must pass an AbortSignal timeout (options.signal).",
          });
        }
      },
    };
  },
};

const noModuleMocks = {
  meta: { name: "no-module-mocks" },
  create(context) {
    return {
      CallExpression(node) {
        const chain = memberChain(node);
        const method = chain.at(-1);
        if (chain[0] === "vi" && (method === "mock" || method === "doMock" || method === "unmock")) {
          context.report({
            node,
            message: `${chain.join(".")} module mocking is banned. Inject collaborators instead.`,
          });
        }
      },
    };
  },
};

const noLeakyMocks = {
  meta: { name: "no-leaky-mocks" },
  create(context) {
    return {
      CallExpression(node) {
        const chain = memberChain(node);
        if (chain[0] === "vi" && (chain.at(-1) === "spyOn" || chain.at(-1) === "fn")) {
          context.report({
            node,
            message:
              "Bare vi.spyOn/vi.fn leaks into the next test. Register cleanup or pass an explicit stub.",
          });
        }
      },
    };
  },
};

const noFakeSuccess = {
  meta: { name: "no-fake-success" },
  create(context) {
    return {
      CallExpression(node) {
        if (!TEST_CALLERS.has(rootCalleeName(node) ?? "")) {
          return;
        }
        const body = node.arguments.at(-1);
        if (
          body &&
          (body.type === "ArrowFunctionExpression" || body.type === "FunctionExpression") &&
          body.body.type === "BlockStatement" &&
          body.body.body.length === 0
        ) {
          context.report({
            node,
            message: "A test with an empty body always passes. Assert observable behavior or delete it.",
          });
        }
      },
    };
  },
};

const noAlwaysTrueAssertion = {
  meta: { name: "no-always-true-assertion" },
  create(context) {
    return {
      CallExpression(node) {
        if (calleeName(node) !== "expect" || node.arguments.length === 0) {
          return;
        }
        const subject = node.arguments[0];
        if (subject.type === "Literal") {
          context.report({
            node,
            message: "Asserting on a literal is always true. Assert on a computed value.",
          });
        }
      },
    };
  },
};

const noDisabledTests = {
  meta: { name: "no-disabled-tests" },
  create(context) {
    return {
      CallExpression(node) {
        const root = rootCalleeName(node);
        if (root === "xit" || root === "xdescribe" || root === "xtest") {
          context.report({ node, message: "Disabled tests are banned. Fix or delete the test." });
          return;
        }
        if (TEST_CALLERS.has(root ?? "") && memberChain(node).includes("skip")) {
          context.report({ node, message: "Disabled tests are banned. Fix or delete the test." });
        }
      },
    };
  },
};

const noFocusedTests = {
  meta: { name: "no-focused-tests" },
  create(context) {
    return {
      CallExpression(node) {
        const root = rootCalleeName(node);
        if (root === "fit" || root === "fdescribe") {
          context.report({ node, message: "Focused tests are banned. They hide the rest of the suite." });
          return;
        }
        if (TEST_CALLERS.has(root ?? "") && memberChain(node).includes("only")) {
          context.report({ node, message: "Focused tests are banned. They hide the rest of the suite." });
        }
      },
    };
  },
};

const preferTestOverIt = {
  meta: { name: "prefer-test-over-it" },
  create(context) {
    return {
      CallExpression(node) {
        if (node.callee.type === "Identifier" && node.callee.name === "it") {
          context.report({ node: node.callee, message: 'Use "test" instead of "it".' });
        }
      },
    };
  },
};

const noUppercaseTestTitle = {
  meta: { name: "no-uppercase-test-title" },
  create(context) {
    return {
      CallExpression(node) {
        if (!TEST_CALLERS.has(rootCalleeName(node) ?? "") || node.arguments.length === 0) {
          return;
        }
        const title = node.arguments[0];
        if (title.type === "Literal" && typeof title.value === "string" && /^[A-Z]/.test(title.value)) {
          context.report({
            node: title,
            message: "Test titles must start lowercase and read as a sentence fragment.",
          });
        }
      },
    };
  },
};

export default {
  meta: { name: "treetime" },
  rules: {
    "no-module-level-mutable": noModuleLevelMutable,
    "no-vague-identifiers": noVagueIdentifiers,
    "no-versioned-names": noVersionedNames,
    "no-typographic-characters": noTypographicCharacters,
    "no-exec-shell-string": noExecShellString,
    "use-class-name-helper": noClassNameHelperBypass,
    "require-io-timeout": requireIoTimeout,
    "no-module-mocks": noModuleMocks,
    "no-leaky-mocks": noLeakyMocks,
    "no-fake-success": noFakeSuccess,
    "no-always-true-assertion": noAlwaysTrueAssertion,
    "no-disabled-tests": noDisabledTests,
    "no-focused-tests": noFocusedTests,
    "prefer-test-over-it": preferTestOverIt,
    "no-uppercase-test-title": noUppercaseTestTitle,
  },
};
