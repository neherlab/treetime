import { defineRule } from "@oxlint/plugins";
import type { ESTree, SourceCode } from "@oxlint/plugins";

import { rootCalleeName, TEST_CALLERS } from "./ast.ts";
import { bindingInitializer, bindingVariable } from "./binding.ts";
import { functionName, isFunctionNode } from "./function-shape.ts";
import type { FunctionNode } from "./function-shape.ts";
import {
  evaluationCompletions,
  observes,
  patternCompletions,
  runtimeExpression,
  staticTruth,
} from "./no-fake-success-evaluation.ts";
import type { ExpressionCompletion, ReturnAnalysis } from "./no-fake-success-evaluation.ts";

const FAILURE_KEYS = new Set(["error", "errors", "cause", "failed", "failure", "reason"]);

const STATUS_KEYS = new Set(["ok", "success", "valid"]);

const PROMISE_FACTORIES = new Set(["resolve", "reject", "all", "allSettled", "any", "race", "try"]);

const PROMISE_CHAINS = new Set(["catch", "then", "finally"]);

export const noFakeSuccessRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description:
        "Disallow functions that ignore their inputs and return a fixed success value, and error handlers that turn an error into success.",
    },
    messages: {
      empty: "A test with an empty body always passes. Assert observable behavior or delete it.",
      ignoredInputs:
        "`{{name}}` never reads its parameters and returns a fixed success value. Implement the path or fail explicitly.",
      swallowedError:
        "This error handler discards the error and reports success. Recover from a specific failure, rethrow with `cause`, or return a failure value.",
    },
  },
  createOnce(context) {
    function checkFunction(node: FunctionNode): void {
      if (node.body === null) {
        return;
      }

      const handler = isErrorHandler(node, context.sourceCode);

      if (handler === undefined) {
        return;
      }

      const parameters = context.sourceCode
        .getDeclaredVariables(node)
        .filter((variable) => variable.defs[0]?.type === "Parameter")
        .filter((variable) => handler || !variable.name.startsWith("_"));

      if (!handler && parameters.length === 0) {
        return;
      }

      const analysis = {
        owner: node,
        parameters,
        sourceCode: context.sourceCode,
        switchEntries: new Map<ESTree.SwitchStatement, number>(),
        staticBlockCompletions,
      };

      if (!handler && observes(node.body, analysis)) {
        return;
      }

      for (const result of successReturns(node.body, analysis)) {
        context.report({
          node: result,
          messageId: handler ? "swallowedError" : "ignoredInputs",
          data: { name: functionName(node, context.sourceCode) ?? "This function" },
        });
      }
    }

    return {
      CallExpression(node) {
        if (!TEST_CALLERS.has(rootCalleeName(node) ?? "")) {
          return;
        }

        const body = node.arguments.at(-1);

        if (
          body !== undefined &&
          (body.type === "ArrowFunctionExpression" || body.type === "FunctionExpression") &&
          body.body !== null &&
          body.body.type === "BlockStatement" &&
          body.body.body.length === 0
        ) {
          context.report({ node, messageId: "empty" });
        }
      },
      ArrowFunctionExpression: checkFunction,
      FunctionDeclaration: checkFunction,
      FunctionExpression: checkFunction,
      CatchClause(node) {
        const analysis = {
          owner: node,
          parameters: context.sourceCode.getDeclaredVariables(node),
          sourceCode: context.sourceCode,
          switchEntries: new Map<ESTree.SwitchStatement, number>(),
          staticBlockCompletions,
        };

        for (const result of successReturns(enclosingBody(node), analysis)) {
          context.report({ node: result, messageId: "swallowedError" });
        }
      },
    };
  },
});

function isErrorHandler(node: FunctionNode, sourceCode: SourceCode): boolean | undefined {
  const { parent } = node;

  if (parent.type !== "CallExpression" || parent.callee.type !== "MemberExpression") {
    return false;
  }

  const { callee } = parent;
  const method = memberName(callee);

  const position =
    (method === "catch" && parent.arguments[0] === node) || (method === "then" && parent.arguments[1] === node);

  return position ? isPromise(callee.object, sourceCode) || undefined : false;
}

function isPromise(
  node: ESTree.Expression,
  sourceCode: SourceCode,
  seen: ReadonlySet<ESTree.Node> = new Set(),
): boolean {
  if (seen.has(node)) {
    return false;
  }

  const visited = new Set([...seen, node]);

  if (node.type === "Identifier") {
    const initializer = bindingInitializer(node, sourceCode);

    return initializer !== undefined && isPromise(initializer, sourceCode, visited);
  }

  if (node.type === "NewExpression") {
    return isPromiseConstructor(node.callee, sourceCode);
  }

  if (node.type !== "CallExpression") {
    return false;
  }

  const { callee } = node;

  if (isFunctionNode(callee)) {
    return callee.async && !callee.generator;
  }

  if (callee.type === "Identifier") {
    if (callee.name === "fetch" && sourceCode.isGlobalReference(callee)) {
      return true;
    }

    const variable = bindingVariable(callee, sourceCode);

    if (variable?.references.some((reference) => reference.isWrite() && !reference.init) === true) {
      return false;
    }

    const definition = variable?.defs[0];

    const fn = definition?.type === "FunctionName" ? definition.node : bindingInitializer(callee, sourceCode);

    return fn !== undefined && isFunctionNode(fn) && fn.async && !fn.generator;
  }

  if (callee.type !== "MemberExpression") {
    return false;
  }

  const method = memberName(callee);

  return (
    method !== undefined &&
    ((PROMISE_FACTORIES.has(method) && isPromiseConstructor(callee.object, sourceCode)) ||
      (PROMISE_CHAINS.has(method) && isPromise(callee.object, sourceCode, visited)))
  );
}

function isPromiseConstructor(node: ESTree.Expression, sourceCode: SourceCode): boolean {
  return node.type === "Identifier" && node.name === "Promise" && sourceCode.isGlobalReference(node);
}

function memberName(node: ESTree.MemberExpression): string | undefined {
  if (!node.computed && node.property.type === "Identifier") {
    return node.property.name;
  }

  return node.property.type === "Literal" && typeof node.property.value === "string" ? node.property.value : undefined;
}

function successReturns(
  body: ESTree.Expression | ESTree.BlockStatement | ESTree.Program,
  analysis: ReturnAnalysis,
): Set<ESTree.Node> {
  const paths =
    body.type === "BlockStatement" || body.type === "Program"
      ? blockCompletions(body.body, false, analysis)
      : expressionCompletions(body, body, false, analysis);

  return new Set(
    paths.flatMap((path) =>
      path.kind === "return" && !path.handled && resultOwner(path.node) === analysis.owner && isSuccessValue(path.value)
        ? [path.node]
        : [],
    ),
  );
}

function enclosingBody(node: ESTree.CatchClause): ESTree.Expression | ESTree.BlockStatement | ESTree.Program {
  let current: ESTree.Node = node.parent;

  while (current.type !== "Program") {
    if (isFunctionNode(current) && current.body !== null) {
      return current.body;
    }

    current = current.parent;
  }

  return current;
}

function resultOwner(node: ESTree.Node): ESTree.Node | undefined {
  let current = node.parent;

  while (current !== null) {
    if (isFunctionNode(current) || current.type === "CatchClause") {
      return current;
    }

    current = current.parent;
  }

  return undefined;
}

type Completion =
  | ExpressionCompletion
  | { kind: "return"; handled: boolean; node: ESTree.Node; value: ESTree.Expression | undefined }
  | { kind: "break" | "continue"; handled: boolean; label: string | undefined };

function statementCompletions(statement: ESTree.Statement, handled: boolean, analysis: ReturnAnalysis): Completion[] {
  if (statement.type === "ReturnStatement") {
    return expressionCompletions(statement.argument ?? undefined, statement, handled, analysis);
  }

  if (statement.type === "ThrowStatement") {
    return evaluationCompletions(statement.argument, handled, analysis, "propagate").map((path) => ({
      kind: "throw",
      handled: path.handled,
    }));
  }

  if (statement.type === "BlockStatement") {
    return blockCompletions(statement.body, handled, analysis);
  }

  if (statement.type === "IfStatement") {
    const truth = staticTruth(statement.test);

    return continueCompletions(evaluationCompletions(statement.test, handled, analysis), (observed) => {
      const alternate: Completion[] =
        statement.alternate === null
          ? [{ kind: "normal", handled: observed } satisfies Completion]
          : statementCompletions(statement.alternate, observed, analysis);

      return [
        ...(truth === false ? [] : statementCompletions(statement.consequent, observed, analysis)),
        ...(truth === true ? [] : alternate),
      ];
    });
  }

  if (statement.type === "TryStatement") {
    return tryCompletions(statement, handled, analysis);
  }

  if (statement.type === "SwitchStatement") {
    return switchCompletions(statement, handled, analysis);
  }

  if (
    statement.type === "ForStatement" ||
    statement.type === "ForInStatement" ||
    statement.type === "ForOfStatement" ||
    statement.type === "WhileStatement" ||
    statement.type === "DoWhileStatement"
  ) {
    return loopCompletions(statement, handled, analysis);
  }

  if (statement.type === "BreakStatement" || statement.type === "ContinueStatement") {
    return [
      {
        kind: statement.type === "BreakStatement" ? "break" : "continue",
        handled,
        label: statement.label?.name,
      },
    ];
  }

  if (statement.type === "LabeledStatement") {
    return statementCompletions(statement.body, handled, analysis).map(
      (path): Completion =>
        path.kind === "break" && path.label === statement.label.name ? { kind: "normal", handled: path.handled } : path,
    );
  }

  if (statement.type === "FunctionDeclaration" || statement.type === "EmptyStatement") {
    return [{ kind: "normal", handled }];
  }

  if (statement.type === "ExpressionStatement") {
    return evaluationCompletions(statement.expression, handled, analysis);
  }

  if (statement.type === "VariableDeclaration") {
    return declarationCompletions(statement, handled, analysis);
  }

  if (statement.type === "ClassDeclaration") {
    return evaluationCompletions(statement, handled, analysis);
  }

  const observed = handled || observes(statement, analysis);

  return [
    { kind: "normal", handled: observed },
    ...(statement.type !== "TSTypeAliasDeclaration" && statement.type !== "TSInterfaceDeclaration"
      ? [{ kind: "throw", handled: observed } satisfies Completion]
      : []),
  ];
}

function continueCompletions(paths: Completion[], next: (handled: boolean) => Completion[]): Completion[] {
  return mergeCompletions(paths.flatMap((path) => (path.kind === "normal" ? next(path.handled) : [path])));
}

function declarationCompletions(
  statement: ESTree.VariableDeclaration,
  handled: boolean,
  analysis: ReturnAnalysis,
): Completion[] {
  let paths: Completion[] = [{ kind: "normal", handled }];

  for (const declaration of statement.declarations) {
    paths = continueCompletions(paths, (observed) =>
      continueCompletions(evaluationCompletions(declaration.init ?? undefined, observed, analysis), (initialized) =>
        patternCompletions(declaration.id, initialized, analysis),
      ),
    );
  }

  return paths;
}

function blockCompletions(statements: ESTree.Statement[], handled: boolean, analysis: ReturnAnalysis): Completion[] {
  let paths: Completion[] = [{ kind: "normal", handled }];

  for (const statement of statements) {
    paths = mergeCompletions(
      paths.flatMap((path) =>
        path.kind === "normal" ? statementCompletions(statement, path.handled, analysis) : [path],
      ),
    );
  }

  return paths;
}

function staticBlockCompletions(
  block: ESTree.StaticBlock,
  handled: boolean,
  analysis: ReturnAnalysis,
): ExpressionCompletion[] {
  return blockCompletions(block.body, handled, analysis).map((path) => {
    if (path.kind !== "normal" && path.kind !== "throw") {
      throw new Error(`Unexpected ${path.kind} completion in a class static block`);
    }

    return path;
  });
}

function mergeCompletions(paths: Completion[]): Completion[] {
  return [...new Map(paths.map((path) => [completionKey(path), path])).values()];
}

function completionKey(path: Completion): string {
  if (path.kind === "return") {
    return `${path.kind}:${path.handled}:${path.node.start}:${path.value?.start}`;
  }

  if (path.kind === "break" || path.kind === "continue") {
    return `${path.kind}:${path.handled}:${path.label ?? ""}`;
  }

  return `${path.kind}:${path.handled}`;
}

function tryCompletions(statement: ESTree.TryStatement, handled: boolean, analysis: ReturnAnalysis): Completion[] {
  const caught = statementCompletions(statement.block, handled, analysis).flatMap((path) =>
    path.kind === "throw" && statement.handler !== null
      ? catchCompletions(statement.handler, path.handled, analysis)
      : [path],
  );

  const finalizer = statement.finalizer;

  if (finalizer === null) {
    return caught;
  }

  return caught.flatMap((path) =>
    statementCompletions(finalizer, path.handled, analysis).map(
      (final): Completion => (final.kind === "normal" ? completionHandled(path, final.handled) : final),
    ),
  );
}

function catchCompletions(handler: ESTree.CatchClause, handled: boolean, analysis: ReturnAnalysis): Completion[] {
  const entry =
    handler.param === null
      ? [{ kind: "normal", handled } satisfies Completion]
      : patternCompletions(handler.param, handled, analysis);

  return continueCompletions(entry, (observed) => statementCompletions(handler.body, observed, analysis));
}

function completionHandled(path: Completion, handled: boolean): Completion {
  if (path.kind === "return") {
    return { kind: path.kind, node: path.node, value: path.value, handled };
  }

  if (path.kind === "break" || path.kind === "continue") {
    return { kind: path.kind, label: path.label, handled };
  }

  return { kind: path.kind, handled };
}

function loopCompletions(statement: LoopStatement, handled: boolean, analysis: ReturnAnalysis): Completion[] {
  const entry = loopEntryCompletions(statement, handled, analysis);

  return continueCompletions(entry, (initialized) => {
    const pending = [initialized];
    const visited = new Set<boolean>();
    const exits: Completion[] = [];

    while (pending.length > 0) {
      const observed = pending.shift();

      if (observed === undefined || visited.has(observed)) {
        continue;
      }

      visited.add(observed);

      const start =
        statement.type === "DoWhileStatement"
          ? [{ kind: "normal", handled: observed } satisfies Completion]
          : loopTestCompletions(statement, observed, analysis);

      for (const path of start) {
        if (path.kind !== "normal") {
          exits.push(path);
          continue;
        }

        const truth = loopTruth(statement);

        if (statement.type !== "DoWhileStatement" && truth !== true) {
          exits.push(path);
        }

        if (statement.type !== "DoWhileStatement" && truth === false) {
          continue;
        }

        const body = continueCompletions(loopBindingCompletions(statement, path.handled, analysis), (bound) =>
          statementCompletions(statement.body, bound, analysis),
        );

        for (const result of body) {
          if (result.kind === "break" && loopOwnsLabel(statement, result.label)) {
            exits.push({ kind: "normal", handled: result.handled });
          } else if (
            result.kind === "normal" ||
            (result.kind === "continue" && loopOwnsLabel(statement, result.label))
          ) {
            const advanced = loopAdvanceCompletions(statement, result.handled, analysis);

            for (const next of advanced) {
              if (next.kind !== "normal") {
                exits.push(next);
              } else {
                if (statement.type === "DoWhileStatement" && truth !== true) {
                  exits.push(next);
                }

                if (statement.type !== "DoWhileStatement" || truth !== false) {
                  pending.push(next.handled);
                }
              }
            }
          } else {
            exits.push(result);
          }
        }
      }
    }

    return mergeCompletions(exits);
  });
}

type LoopStatement =
  | ESTree.ForStatement
  | ESTree.ForInStatement
  | ESTree.ForOfStatement
  | ESTree.WhileStatement
  | ESTree.DoWhileStatement;

function loopBindingCompletions(statement: LoopStatement, handled: boolean, analysis: ReturnAnalysis): Completion[] {
  if (statement.type !== "ForInStatement" && statement.type !== "ForOfStatement") {
    return [{ kind: "normal", handled }];
  }

  if (statement.left.type !== "VariableDeclaration") {
    return patternCompletions(statement.left, handled, analysis);
  }

  const declaration = statement.left.declarations[0];

  return declaration === undefined
    ? [{ kind: "normal", handled }]
    : patternCompletions(declaration.id, handled, analysis);
}

function loopAdvanceCompletions(statement: LoopStatement, handled: boolean, analysis: ReturnAnalysis): Completion[] {
  if (statement.type === "ForStatement") {
    return evaluationCompletions(statement.update ?? undefined, handled, analysis);
  }

  return statement.type === "DoWhileStatement"
    ? loopTestCompletions(statement, handled, analysis)
    : [{ kind: "normal", handled }];
}

function loopEntryCompletions(statement: LoopStatement, handled: boolean, analysis: ReturnAnalysis): Completion[] {
  if (statement.type === "ForInStatement" || statement.type === "ForOfStatement") {
    return evaluationCompletions(statement.right, handled, analysis);
  }

  if (statement.type === "ForStatement" && statement.init !== null) {
    return statement.init.type === "VariableDeclaration"
      ? declarationCompletions(statement.init, handled, analysis)
      : evaluationCompletions(statement.init, handled, analysis);
  }

  return [{ kind: "normal", handled }];
}

function loopTestCompletions(statement: LoopStatement, handled: boolean, analysis: ReturnAnalysis): Completion[] {
  if (statement.type === "ForInStatement" || statement.type === "ForOfStatement") {
    return [
      { kind: "normal", handled },
      ...(loopTruth(statement) === false ? [] : [{ kind: "throw", handled } satisfies Completion]),
    ];
  }

  return evaluationCompletions(statement.test ?? undefined, handled, analysis);
}

function loopTruth(statement: LoopStatement): boolean | undefined {
  if (
    statement.type === "ForInStatement" &&
    statement.right.type === "Literal" &&
    !("regex" in statement.right) &&
    statement.right.value === null
  ) {
    return false;
  }

  if (statement.type === "ForInStatement" || statement.type === "ForOfStatement") {
    return undefined;
  }

  return statement.test === null ? true : staticTruth(statement.test);
}

function loopOwnsLabel(statement: LoopStatement, label: string | undefined): boolean {
  if (label === undefined) {
    return true;
  }

  let parent = statement.parent;

  while (parent.type === "LabeledStatement") {
    if (parent.label.name === label) {
      return true;
    }

    parent = parent.parent;
  }

  return false;
}

function switchCompletions(
  statement: ESTree.SwitchStatement,
  handled: boolean,
  analysis: ReturnAnalysis,
): Completion[] {
  let unmatched: Completion[] = evaluationCompletions(statement.discriminant, handled, analysis);
  const branches: Completion[] = [];

  for (const [index, branch] of statement.cases.entries()) {
    const test = branch.test;

    if (test === null) {
      continue;
    }

    unmatched = continueCompletions(unmatched, (observed) => evaluationCompletions(test, observed, analysis));

    const matches =
      statement.discriminant.type === "Literal" &&
      test.type === "Literal" &&
      !("regex" in statement.discriminant) &&
      !("regex" in test)
        ? statement.discriminant.value === test.value
        : undefined;

    if (matches !== false) {
      branches.push(
        ...continueCompletions(unmatched, (observed) => switchBranchCompletions(statement, index, observed, analysis)),
      );
    }

    if (matches === true) {
      unmatched = unmatched.filter((path) => path.kind !== "normal");
    }
  }

  const fallback = statement.cases.findIndex((branch) => branch.test === null);
  branches.push(
    ...continueCompletions(unmatched, (observed) =>
      fallback < 0
        ? [{ kind: "normal", handled: observed }]
        : switchBranchCompletions(statement, fallback, observed, analysis),
    ),
  );

  return mergeCompletions(branches);
}

function switchBranchCompletions(
  statement: ESTree.SwitchStatement,
  index: number,
  handled: boolean,
  analysis: ReturnAnalysis,
): Completion[] {
  return blockCompletions(
    statement.cases.slice(index).flatMap((branch) => branch.consequent),
    handled,
    { ...analysis, switchEntries: new Map([...analysis.switchEntries, [statement, index]]) },
  ).map(
    (path): Completion =>
      path.kind === "break" && path.label === undefined ? { kind: "normal", handled: path.handled } : path,
  );
}

function expressionCompletions(
  expression: ESTree.Expression | undefined,
  node: ESTree.Node,
  handled: boolean,
  analysis: ReturnAnalysis,
): Completion[] {
  if (expression?.type === "ConditionalExpression") {
    const truth = staticTruth(expression.test);

    return continueCompletions(evaluationCompletions(expression.test, handled, analysis), (observed) => [
      ...(truth === false ? [] : expressionCompletions(expression.consequent, node, observed, analysis)),
      ...(truth === true ? [] : expressionCompletions(expression.alternate, node, observed, analysis)),
    ]);
  }

  if (expression?.type === "SequenceExpression") {
    let paths: Completion[] = [{ kind: "normal", handled }];

    for (const part of expression.expressions.slice(0, -1)) {
      paths = continueCompletions(paths, (observed) => evaluationCompletions(part, observed, analysis));
    }

    return continueCompletions(paths, (observed) =>
      expressionCompletions(expression.expressions.at(-1), node, observed, analysis),
    );
  }

  const runtime = expression === undefined ? undefined : runtimeExpression(expression);

  if (runtime !== expression) {
    return expressionCompletions(runtime, node, handled, analysis);
  }

  return evaluationCompletions(expression, handled, analysis, "propagate").map(
    (path): Completion =>
      path.kind === "normal" ? { kind: "return", node, value: expression, handled: path.handled } : path,
  );
}

function isSuccessValue(expression: ESTree.Expression | undefined): boolean {
  if (expression === undefined) {
    return false;
  }

  if (expression.type === "Literal") {
    return expression.value === true;
  }

  if (expression.type === "ArrayExpression") {
    return expression.elements.every(isPlainElement);
  }

  if (expression.type === "ObjectExpression") {
    return expression.properties.some(isSuccessStatus) && expression.properties.every(isSuccessProperty);
  }

  return false;
}

function isSuccessStatus(property: ESTree.ObjectPropertyKind): boolean {
  if (property.type !== "Property") {
    return false;
  }

  const name = propertyName(property);

  return (
    name !== undefined && STATUS_KEYS.has(name) && property.value.type === "Literal" && property.value.value === true
  );
}

function isSuccessProperty(property: ESTree.ObjectPropertyKind): boolean {
  if (property.type !== "Property") {
    return false;
  }

  const name = propertyName(property);

  if (name === undefined) {
    return false;
  }

  if (FAILURE_KEYS.has(name)) {
    return false;
  }

  if (STATUS_KEYS.has(name)) {
    return property.value.type === "Literal" && property.value.value === true;
  }

  return isPlainValue(property.value);
}

function propertyName(property: ESTree.ObjectProperty): string | undefined {
  if (property.key.type === "Literal" && typeof property.key.value === "string") {
    return property.key.value;
  }

  return !property.computed && property.key.type === "Identifier" ? property.key.name : undefined;
}

function isPlainElement(element: ESTree.ArrayExpressionElement): boolean {
  return element !== null && element.type !== "SpreadElement" && isPlainValue(element);
}

function isPlainValue(node: ESTree.Expression): boolean {
  if (node.type === "Literal") {
    return true;
  }

  if (node.type === "TemplateLiteral") {
    return node.expressions.length === 0;
  }

  if (node.type === "UnaryExpression") {
    return node.argument.type === "Literal";
  }

  if (node.type === "ArrayExpression") {
    return node.elements.every(isPlainElement);
  }

  if (node.type === "ObjectExpression") {
    return node.properties.every(isSuccessProperty);
  }

  return false;
}
