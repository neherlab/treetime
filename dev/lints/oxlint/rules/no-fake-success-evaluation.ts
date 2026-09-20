import type { ESTree, SourceCode, Variable } from "@oxlint/plugins";

import { bindingVariable } from "./binding.ts";
import { enclosingFunction, isFunctionNode } from "./function-shape.ts";
import type { FunctionNode } from "./function-shape.ts";

export function evaluationCompletions(
  expression: ESTree.Expression | undefined,
  handled: boolean,
  analysis: ReturnAnalysis,
  usage: ValueUsage = "observe",
): ExpressionCompletion[] {
  return evaluate(expression, handled, analysis, usage).map((path) => ({
    kind: path.kind === "throw" ? "throw" : "normal",
    handled: path.handled,
  }));
}

export function patternCompletions(
  pattern: ESTree.Node,
  handled: boolean,
  analysis: ReturnAnalysis,
): ExpressionCompletion[] {
  return evaluatePattern(pattern, handled, analysis).map((path) => ({
    kind: path.kind === "throw" ? "throw" : "normal",
    handled: path.handled,
  }));
}

export function staticTruth(expression: ESTree.Expression): boolean | undefined {
  const runtime = expressionResult(expression);
  if (runtime !== expression) {
    return staticTruth(runtime);
  }
  if (expression.type === "Literal") {
    return "regex" in expression || Boolean(expression.value);
  }
  if (expression.type === "UnaryExpression" && expression.operator === "!") {
    const truth = staticTruth(expression.argument);
    return truth === undefined ? undefined : !truth;
  }
  if (expression.type === "LogicalExpression") {
    const right = logicalRightEvaluated(expression);
    if (right === undefined) {
      return undefined;
    }
    return staticTruth(right ? expression.right : expression.left);
  }
  return undefined;
}

export function runtimeExpression(expression: ESTree.Expression): ESTree.Expression {
  if (
    expression.type === "TSAsExpression" ||
    expression.type === "TSSatisfiesExpression" ||
    expression.type === "TSNonNullExpression" ||
    expression.type === "TSTypeAssertion" ||
    expression.type === "TSInstantiationExpression" ||
    expression.type === "ParenthesizedExpression"
  ) {
    return runtimeExpression(expression.expression);
  }
  return expression;
}

export function observes(region: ESTree.Node, analysis: ReturnAnalysis): boolean {
  return analysis.parameters.some((variable) =>
    variable.references.some(
      (reference) => reference.isRead() && isActiveReference(reference.identifier, region, analysis),
    ),
  );
}

export interface ReturnAnalysis {
  owner: ESTree.Node;
  parameters: Variable[];
  sourceCode: SourceCode;
  switchEntries: ReadonlyMap<ESTree.SwitchStatement, number>;
  staticBlockCompletions: (
    block: ESTree.StaticBlock,
    handled: boolean,
    analysis: ReturnAnalysis,
  ) => ExpressionCompletion[];
}

export interface ExpressionCompletion {
  kind: "normal" | "throw";
  handled: boolean;
}

type ValueUsage = "observe" | "propagate" | "callee";

interface Evaluation {
  kind: "normal" | "throw" | "shortCircuit";
  handled: boolean;
  callable: boolean;
}

function evaluate(
  expression: ESTree.Expression | undefined,
  handled: boolean,
  analysis: ReturnAnalysis,
  usage: ValueUsage = "observe",
): Evaluation[] {
  if (expression === undefined) {
    return normal(handled);
  }
  if (isFunctionNode(expression)) {
    return [
      {
        kind: "normal",
        handled,
        callable: expression.body !== null && observes(expression.body, analysis),
      },
    ];
  }
  if (expression.type === "Identifier") {
    return evaluateIdentifier(expression, handled, analysis, usage);
  }
  const runtime = runtimeExpression(expression);
  if (runtime !== expression) {
    return evaluate(runtime, handled, analysis, usage);
  }
  if (expression.type === "ChainExpression") {
    return evaluate(expression.expression, handled, analysis, usage).map((path) =>
      path.kind === "shortCircuit" ? { kind: "normal", handled: path.handled, callable: path.callable } : path,
    );
  }
  if (expression.type === "ConditionalExpression") {
    return evaluateConditional(expression, handled, analysis, usage);
  }
  if (expression.type === "LogicalExpression") {
    return evaluateLogical(expression, handled, analysis, usage);
  }
  if (expression.type === "SequenceExpression") {
    return after(evaluateList(expression.expressions.slice(0, -1), handled, analysis), (path) =>
      evaluate(expression.expressions.at(-1), path.handled, analysis, usage),
    );
  }
  if (expression.type === "CallExpression" || expression.type === "NewExpression") {
    return evaluateCall(expression, handled, analysis);
  }
  if (expression.type === "MemberExpression") {
    return after(evaluateMemberReference(expression, handled, analysis), (path) => mayThrow(path.handled));
  }
  if (expression.type === "ArrayExpression") {
    return evaluateList(expression.elements, handled, analysis, usage);
  }
  if (expression.type === "ObjectExpression") {
    return evaluateObject(expression, handled, analysis, usage);
  }
  if (expression.type === "TemplateLiteral") {
    return evaluateTemplate(expression, handled, analysis);
  }
  if (expression.type === "TaggedTemplateExpression") {
    return after(evaluate(expression.tag, handled, analysis, "callee"), (tag) =>
      after(evaluateArguments(expression.quasi.expressions, tag, analysis), (path) =>
        mayThrow(path.handled || path.callable),
      ),
    );
  }
  if (expression.type === "AssignmentExpression") {
    return evaluateAssignment(expression, handled, analysis);
  }
  if (expression.type === "ClassExpression" || expression.type === "ClassDeclaration") {
    return evaluateClass(expression, handled, analysis);
  }
  return evaluateOperation(expression, handled, analysis);
}

function evaluateClass(expression: ESTree.Class, handled: boolean, analysis: ReturnAnalysis): Evaluation[] {
  if (expression.declare === true) {
    return normal(handled);
  }
  let paths =
    expression.superClass === null
      ? normal(handled)
      : after(evaluate(expression.superClass, handled, analysis), (path) =>
          classHeritageCheck(expression.superClass, path.handled),
        );
  for (const element of expression.body.body) {
    if (
      element.type !== "StaticBlock" &&
      element.type !== "TSIndexSignature" &&
      element.computed &&
      element.key.type !== "PrivateIdentifier"
    ) {
      const key = element.key;
      paths = after(paths, (path) => evaluateKey(key, path.handled, analysis));
    }
    if (element.type === "MethodDefinition" && element.static) {
      paths = after(paths, (path) => classStaticPropertyDefinition(element, path.handled));
    }
  }
  for (const element of expression.body.body) {
    paths = after(paths, (path) => {
      if (element.type === "StaticBlock") {
        return analysis
          .staticBlockCompletions(element, path.handled, analysis)
          .map((result) => ({ kind: result.kind, handled: result.handled, callable: false }));
      }
      if (
        (element.type === "PropertyDefinition" || element.type === "AccessorProperty") &&
        element.static &&
        element.declare !== true
      ) {
        return after(evaluate(element.value ?? undefined, path.handled, analysis), (initialized) =>
          classStaticPropertyDefinition(element, initialized.handled),
        );
      }
      return normal(path.handled);
    });
  }
  return paths.map((path) => ({ kind: path.kind, handled: path.handled, callable: false }));
}

function classHeritageCheck(expression: ESTree.Expression | null, handled: boolean): Evaluation[] {
  if (expression === null) {
    return normal(handled);
  }
  const runtime = expressionResult(expression);
  if (runtime.type === "Literal") {
    return !("regex" in runtime) && runtime.value === null
      ? normal(handled)
      : [{ kind: "throw", handled, callable: false }];
  }
  if (
    runtime.type === "ClassExpression" ||
    (runtime.type === "FunctionExpression" && !runtime.async && !runtime.generator)
  ) {
    return normal(handled);
  }
  if (runtime.type === "ArrowFunctionExpression" || runtime.type === "FunctionExpression") {
    return [{ kind: "throw", handled, callable: false }];
  }
  return mayThrow(handled);
}

function classStaticPropertyDefinition(
  element: ESTree.MethodDefinition | ESTree.PropertyDefinition | ESTree.AccessorProperty,
  handled: boolean,
): Evaluation[] {
  const key = element.key.type === "PrivateIdentifier" ? element.key : expressionResult(element.key);
  if (key.type === "Literal") {
    return key.value === "prototype" ? [{ kind: "throw", handled, callable: false }] : normal(handled);
  }
  return element.computed ? mayThrow(handled) : normal(handled);
}

function evaluateIdentifier(
  expression: Extract<ESTree.Expression, { type: "Identifier" }>,
  handled: boolean,
  analysis: ReturnAnalysis,
  usage: ValueUsage,
): Evaluation[] {
  const read = analysis.parameters.some((variable) =>
    variable.references.some((reference) => reference.isRead() && reference.identifier === expression),
  );
  const initialized = identifierInitialized(expression, analysis, usage);
  if (initialized === false) {
    return [{ kind: "throw", handled, callable: false }];
  }
  return [
    {
      kind: "normal",
      handled: handled || (usage !== "propagate" && read),
      callable: !read && observes(expression, analysis),
    },
    ...(initialized === undefined ? [{ kind: "throw", handled, callable: false } satisfies Evaluation] : []),
  ];
}

function evaluateConditional(
  expression: ESTree.ConditionalExpression,
  handled: boolean,
  analysis: ReturnAnalysis,
  usage: ValueUsage,
): Evaluation[] {
  const truth = staticTruth(expression.test);
  return after(evaluate(expression.test, handled, analysis), (path) => [
    ...(truth === false ? [] : evaluate(expression.consequent, path.handled, analysis, usage)),
    ...(truth === true ? [] : evaluate(expression.alternate, path.handled, analysis, usage)),
  ]);
}

function evaluateLogical(
  expression: ESTree.LogicalExpression,
  handled: boolean,
  analysis: ReturnAnalysis,
  usage: ValueUsage,
): Evaluation[] {
  const right = logicalRightEvaluated(expression);
  return after(evaluate(expression.left, handled, analysis, usage), (path) => [
    ...(right === true ? [] : [path]),
    ...(right === false ? [] : evaluate(expression.right, path.handled, analysis, usage)),
  ]);
}

function evaluateCall(
  expression: ESTree.CallExpression | ESTree.NewExpression,
  handled: boolean,
  analysis: ReturnAnalysis,
): Evaluation[] {
  const callee = evaluate(expression.callee, handled, analysis, "callee");
  const selected =
    expression.type === "CallExpression" && expression.optional ? optionalPaths(callee, expression.callee) : callee;
  return after(selected, (path) =>
    after(evaluateArguments(expression.arguments, path, analysis), (argumentsPath) =>
      mayThrow(argumentsPath.handled || argumentsPath.callable),
    ),
  );
}

function evaluateArguments(
  expressions: (ESTree.Expression | ESTree.SpreadElement)[],
  callee: Evaluation,
  analysis: ReturnAnalysis,
): Evaluation[] {
  let paths = [callee];
  for (const expression of expressions) {
    paths = after(paths, (before) =>
      evaluateElement(expression, before.handled, analysis, "observe").map((path) => ({
        kind: path.kind,
        handled: path.handled,
        callable: before.callable || path.callable,
      })),
    );
  }
  return paths;
}

function evaluateMemberReference(
  expression: ESTree.MemberExpression,
  handled: boolean,
  analysis: ReturnAnalysis,
): Evaluation[] {
  const object = evaluate(expression.object, handled, analysis);
  const selected = expression.optional ? optionalPaths(object, expression.object) : object;
  return after(selected, (path) =>
    expression.computed ? evaluate(expression.property, path.handled, analysis) : normal(path.handled),
  );
}

function optionalPaths(paths: Evaluation[], base: ESTree.Expression): Evaluation[] {
  const nullish = staticNullish(base);
  return after(paths, (path) => [
    ...(nullish === true ? [] : [path]),
    ...(nullish === false
      ? []
      : [{ kind: "shortCircuit", handled: path.handled, callable: false } satisfies Evaluation]),
  ]);
}

function evaluateList(
  expressions: (ESTree.Expression | ESTree.SpreadElement | null)[],
  handled: boolean,
  analysis: ReturnAnalysis,
  usage: ValueUsage = "observe",
): Evaluation[] {
  let paths = normal(handled);
  for (const expression of expressions) {
    paths = after(paths, (path) => evaluateElement(expression, path.handled, analysis, usage));
  }
  return paths.map((path) => ({ kind: path.kind, handled: path.handled, callable: false }));
}

function evaluateElement(
  expression: ESTree.Expression | ESTree.SpreadElement | null,
  handled: boolean,
  analysis: ReturnAnalysis,
  usage: ValueUsage,
): Evaluation[] {
  if (expression === null) {
    return normal(handled);
  }
  if (expression.type === "SpreadElement") {
    return after(evaluate(expression.argument, handled, analysis, usage), (path) => mayThrow(path.handled));
  }
  return evaluate(expression, handled, analysis, usage);
}

function evaluateObject(
  expression: ESTree.ObjectExpression,
  handled: boolean,
  analysis: ReturnAnalysis,
  usage: ValueUsage,
): Evaluation[] {
  let paths = normal(handled);
  for (const property of expression.properties) {
    paths = after(paths, (before) => {
      if (property.type === "SpreadElement") {
        return evaluateElement(property, before.handled, analysis, usage);
      }
      const key =
        property.computed && property.key.type !== "PrivateIdentifier"
          ? evaluateKey(property.key, before.handled, analysis)
          : normal(before.handled);
      return after(key, (path) => evaluate(property.value, path.handled, analysis, usage));
    });
  }
  return paths.map((path) => ({ kind: path.kind, handled: path.handled, callable: false }));
}

function evaluateKey(key: ESTree.Expression, handled: boolean, analysis: ReturnAnalysis): Evaluation[] {
  const result = expressionResult(key);
  return after(evaluate(key, handled, analysis), (path) =>
    result.type === "Literal" && !("regex" in result) ? normal(path.handled) : mayThrow(path.handled),
  );
}

function evaluateTemplate(
  expression: ESTree.TemplateLiteral,
  handled: boolean,
  analysis: ReturnAnalysis,
): Evaluation[] {
  let paths = normal(handled);
  for (const substitution of expression.expressions) {
    paths = after(paths, (before) =>
      after(evaluate(substitution, before.handled, analysis), (path) =>
        substitution.type === "Literal" && !("regex" in substitution) ? normal(path.handled) : mayThrow(path.handled),
      ),
    );
  }
  return paths;
}

function evaluateAssignment(
  expression: ESTree.AssignmentExpression,
  handled: boolean,
  analysis: ReturnAnalysis,
): Evaluation[] {
  if (expression.left.type === "ArrayPattern" || expression.left.type === "ObjectPattern") {
    return after(evaluate(expression.right, handled, analysis), (path) =>
      evaluatePattern(expression.left, path.handled, analysis),
    );
  }
  const target = expression.left;
  const reference =
    target.type === "MemberExpression" ? evaluateMemberReference(target, handled, analysis) : normal(handled);
  let left = reference;
  if (expression.operator !== "=") {
    left =
      target.type === "MemberExpression"
        ? after(reference, (path) => mayThrow(path.handled))
        : evaluate(target, handled, analysis);
  }
  return after(left, (before) => [
    ...(["&&=", "||=", "??="].includes(expression.operator) ? normal(before.handled) : []),
    ...after(evaluate(expression.right, before.handled, analysis), (path) => mayThrow(path.handled)),
  ]);
}

function evaluatePattern(pattern: ESTree.Node, handled: boolean, analysis: ReturnAnalysis): Evaluation[] {
  if (pattern.type === "Identifier") {
    return normal(handled);
  }
  if (pattern.type === "RestElement") {
    return evaluatePatternElement(pattern.argument, handled, analysis);
  }
  if (pattern.type === "AssignmentPattern") {
    return after([...normal(handled), ...evaluate(pattern.right, handled, analysis)], (path) =>
      evaluatePattern(pattern.left, path.handled, analysis),
    );
  }
  if (pattern.type === "ObjectPattern") {
    let paths = mayThrow(handled);
    for (const property of pattern.properties) {
      paths = after(paths, (before) => {
        if (property.type === "RestElement") {
          return evaluatePattern(property, before.handled, analysis);
        }
        const key =
          property.computed && property.key.type !== "PrivateIdentifier"
            ? evaluateKey(property.key, before.handled, analysis)
            : normal(before.handled);
        return after(key, (path) => evaluatePatternElement(property.value, path.handled, analysis));
      });
    }
    return paths;
  }
  if (pattern.type === "ArrayPattern") {
    let paths = mayThrow(handled);
    for (const element of pattern.elements) {
      paths = after(paths, (path) =>
        element === null ? mayThrow(path.handled) : evaluatePatternElement(element, path.handled, analysis),
      );
    }
    return after(paths, (path) => mayThrow(path.handled));
  }
  if (pattern.type === "MemberExpression") {
    return after(evaluateMemberReference(pattern, handled, analysis), (path) => mayThrow(path.handled));
  }
  if (
    pattern.type === "TSAsExpression" ||
    pattern.type === "TSSatisfiesExpression" ||
    pattern.type === "TSNonNullExpression" ||
    pattern.type === "TSTypeAssertion"
  ) {
    return evaluatePattern(pattern.expression, handled, analysis);
  }
  return mayThrow(handled);
}

function evaluatePatternElement(element: ESTree.Node, handled: boolean, analysis: ReturnAnalysis): Evaluation[] {
  if (element.type === "RestElement") {
    return evaluatePatternElement(element.argument, handled, analysis);
  }
  const target = element.type === "AssignmentPattern" ? element.left : element;
  const reference =
    target.type === "MemberExpression" ? evaluateMemberReference(target, handled, analysis) : normal(handled);
  return after(reference, (before) =>
    after(mayThrow(before.handled), (value) => {
      const initialized =
        element.type === "AssignmentPattern"
          ? [...normal(value.handled), ...evaluate(element.right, value.handled, analysis)]
          : normal(value.handled);
      return after(initialized, (path) =>
        target.type === "MemberExpression" ? mayThrow(path.handled) : evaluatePattern(target, path.handled, analysis),
      );
    }),
  );
}

function evaluateOperation(expression: ESTree.Expression, handled: boolean, analysis: ReturnAnalysis): Evaluation[] {
  if (expression.type === "BinaryExpression") {
    const left =
      expression.left.type === "PrivateIdentifier" ? normal(handled) : evaluate(expression.left, handled, analysis);
    return after(left, (before) =>
      after(evaluate(expression.right, before.handled, analysis), (path) =>
        expression.operator === "===" || expression.operator === "!==" ? normal(path.handled) : mayThrow(path.handled),
      ),
    );
  }
  if (expression.type === "UnaryExpression") {
    if (
      expression.operator === "typeof" &&
      expression.argument.type === "Identifier" &&
      bindingVariable(expression.argument, analysis.sourceCode) === undefined
    ) {
      return normal(handled);
    }
    return after(evaluate(expression.argument, handled, analysis), (path) =>
      ["!", "void", "typeof"].includes(expression.operator) ? normal(path.handled) : mayThrow(path.handled),
    );
  }
  if (expression.type === "AwaitExpression" || expression.type === "YieldExpression") {
    return after(evaluate(expression.argument ?? undefined, handled, analysis), (path) => mayThrow(path.handled));
  }
  if (expression.type === "UpdateExpression") {
    return after(evaluate(expression.argument, handled, analysis), (path) => mayThrow(path.handled));
  }
  if (expression.type === "ImportExpression") {
    return after(evaluate(expression.source, handled, analysis), (before) =>
      after(evaluate(expression.options ?? undefined, before.handled, analysis), (path) => mayThrow(path.handled)),
    );
  }
  if (expression.type === "Literal" || expression.type === "ThisExpression" || expression.type === "MetaProperty") {
    return normal(handled);
  }
  return mayThrow(handled);
}

function after(paths: Evaluation[], next: (path: Evaluation) => Evaluation[]): Evaluation[] {
  const continued = paths.flatMap((path) => (path.kind === "normal" ? next(path) : [path]));
  return [...new Map(continued.map((path) => [`${path.kind}:${path.handled}:${path.callable}`, path])).values()];
}

function normal(handled: boolean): Evaluation[] {
  return [{ kind: "normal", handled, callable: false }];
}

function mayThrow(handled: boolean): Evaluation[] {
  return [...normal(handled), { kind: "throw", handled, callable: false }];
}

function logicalRightEvaluated(expression: ESTree.LogicalExpression): boolean | undefined {
  if (expression.operator === "??") {
    return staticNullish(expression.left);
  }
  const truth = staticTruth(expression.left);
  if (truth === undefined) {
    return undefined;
  }
  return expression.operator === "&&" ? truth : !truth;
}

function staticNullish(expression: ESTree.Expression): boolean | undefined {
  const runtime = expressionResult(expression);
  if (runtime !== expression) {
    return staticNullish(runtime);
  }
  if (expression.type === "Literal") {
    return !("regex" in expression) && expression.value === null;
  }
  if (expression.type === "UnaryExpression" && expression.operator === "void") {
    return true;
  }
  return undefined;
}

function expressionResult(expression: ESTree.Expression): ESTree.Expression {
  const runtime = runtimeExpression(expression);
  if (runtime.type === "SequenceExpression") {
    const last = runtime.expressions.at(-1);
    return last === undefined ? runtime : expressionResult(last);
  }
  return runtime;
}

function identifierInitialized(
  expression: Extract<ESTree.Expression, { type: "Identifier" }>,
  analysis: ReturnAnalysis,
  usage: ValueUsage,
): boolean | undefined {
  const variable = bindingVariable(expression, analysis.sourceCode);
  const definition = variable?.defs[0];
  if (usage === "callee" && (definition === undefined || definition.type === "ImportBinding")) {
    return true;
  }
  if (variable === undefined || definition === undefined) {
    return undefined;
  }
  if (
    definition.type === "Variable" &&
    definition.parent?.type === "VariableDeclaration" &&
    definition.parent.kind === "var" &&
    (variable.scope.variableScope.type === "function" || variable.scope.variableScope.type === "class-static-block")
  ) {
    return true;
  }
  if (definition.type === "ClassName" && variable.scope.type === "class") {
    return classBindingInitialized(expression, variable.scope.block);
  }
  if (
    definition.type === "FunctionName" ||
    (definition.type === "CatchClause" &&
      definition.node.type === "CatchClause" &&
      isWithin(expression, definition.node.body)) ||
    (definition.type === "Parameter" &&
      isFunctionNode(definition.node) &&
      definition.node.body !== null &&
      isWithin(expression, definition.node.body))
  ) {
    return true;
  }
  if (definition.type !== "Variable" && definition.type !== "ClassName") {
    return undefined;
  }
  if (definition.node.end < expression.start) {
    return declarationExecuted(variable, definition.node, analysis);
  }
  if (definition.node.type === "VariableDeclarator" && isWithin(expression, definition.node.id)) {
    return undefined;
  }
  return enclosingFunction(expression) === enclosingFunction(definition.node) ? false : undefined;
}

function declarationExecuted(
  variable: Variable,
  declaration: ESTree.Node,
  analysis: ReturnAnalysis,
): boolean | undefined {
  const block = variable.scope.block;
  if (block.type !== "SwitchStatement") {
    return true;
  }
  const entry = analysis.switchEntries.get(block);
  if (entry === undefined) {
    return isWithin(block, analysis.owner) ? false : undefined;
  }
  return block.cases.slice(entry).some((branch) => isWithin(declaration, branch));
}

function classBindingInitialized(node: ESTree.Node, owner: ESTree.Node): boolean | undefined {
  let current: ESTree.Node | null = node;
  let deferred = false;
  while (current !== null && current !== owner) {
    const parent: ESTree.Node | null = current.parent;
    if (parent?.type === "MethodDefinition" && parent.value === current && parent.parent.parent === owner) {
      return true;
    }
    if (parent?.type === "ClassBody" && parent.parent === owner && current.type === "StaticBlock") {
      return true;
    }
    if (
      (parent?.type === "PropertyDefinition" || parent?.type === "AccessorProperty") &&
      parent.value === current &&
      parent.parent.parent === owner
    ) {
      return true;
    }
    deferred ||= isFunctionNode(current);
    current = parent;
  }
  return deferred ? undefined : false;
}

function isActiveReference(
  node: ESTree.Node,
  region: ESTree.Node,
  analysis: ReturnAnalysis,
  seen: ReadonlySet<FunctionNode> = new Set(),
): boolean {
  let current: ESTree.Node | null = node;
  while (current !== null && current !== analysis.owner) {
    const parent: ESTree.Node | null = current.parent;
    if (
      (parent?.type === "PropertyDefinition" || parent?.type === "AccessorProperty") &&
      !parent.static &&
      parent.value === current
    ) {
      return false;
    }
    if (isFunctionNode(current)) {
      return functionMayRun(current, region, analysis, seen);
    }
    if (current === region) {
      return true;
    }
    current = current.parent;
  }
  return false;
}

function functionMayRun(
  node: FunctionNode,
  region: ESTree.Node,
  analysis: ReturnAnalysis,
  seen: ReadonlySet<FunctionNode>,
): boolean {
  if (seen.has(node)) {
    return false;
  }
  if (classFunctionDeferred(node)) {
    return false;
  }
  const visited = new Set([...seen, node]);
  const declaration = node.parent.type === "VariableDeclarator" ? node.parent : node;
  if (declaration === node && node.type !== "FunctionDeclaration") {
    return isActiveReference(node.parent, region, analysis, visited);
  }
  return analysis.sourceCode
    .getDeclaredVariables(declaration)
    .some((variable) =>
      variable.references.some(
        (reference) =>
          reference.isRead() &&
          !classFunctionDeferred(reference.identifier) &&
          !isWithin(reference.identifier, node) &&
          isActiveReference(reference.identifier, region, analysis, visited),
      ),
    );
}

function classFunctionDeferred(node: ESTree.Node): boolean {
  let current: ESTree.Node | null = node.parent;
  while (current !== null && !isFunctionNode(current)) {
    if (
      current.type === "MethodDefinition" ||
      current.type === "PropertyDefinition" ||
      current.type === "AccessorProperty" ||
      current.type === "ClassDeclaration" ||
      current.type === "ClassExpression"
    ) {
      return true;
    }
    if (
      current.type === "CallExpression" ||
      current.type === "NewExpression" ||
      current.type === "TaggedTemplateExpression" ||
      current.type === "BlockStatement" ||
      current.type === "StaticBlock"
    ) {
      return false;
    }
    current = current.parent;
  }
  return false;
}

function isWithin(node: ESTree.Node, ancestor: ESTree.Node): boolean {
  let current: ESTree.Node | null = node;
  while (current !== null && current !== ancestor) {
    current = current.parent;
  }
  return current === ancestor;
}
