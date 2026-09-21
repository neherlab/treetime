import { noSideEffectsInGettersRule } from "../rules/no-side-effects-in-getters.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

const error = { messageId: "sideEffect" };

tester.run("treetime/no-side-effects-in-getters", noSideEffectsInGettersRule, {
  valid: [
    "class C { get x() { return this._x } }",
    "class C { get x() { const a = []; a.push(1); return a } }",
    "class C { get handler() { return () => { this.count++ } } }",
    "const o = { get y() { return 1 } }",
    "class C { method() { this._x = 1 } }",
  ],
  invalid: [
    { code: "class C { get x() { this._x = 1; return this._x } }", errors: [error] },
    { code: "class C { get n() { this._n++; return this._n } }", errors: [error] },
    {
      code: "class C { get items() { this._items.push(1); return this._items } }",
      errors: [error],
    },
    {
      code: "let outer = []; const o = { get z() { outer.push(1); return outer } }",
      errors: [error],
    },
  ],
});
