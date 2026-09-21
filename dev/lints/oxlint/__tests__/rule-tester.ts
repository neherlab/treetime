import { describe, it } from "node:test";

import { RuleTester } from "oxlint/plugins-dev";

RuleTester.describe = (name, fn): void => {
  void describe(name, fn);
};

RuleTester.it = (name, fn): void => {
  // oxlint-disable-next-line treetime/prefer-test-over-it -- RuleTester binds the runner's it
  void it(name, fn);
};

export function ruleTester(lang: "ts" | "tsx", globals: Record<string, "readonly"> = {}): RuleTester {
  return new RuleTester({ languageOptions: { parserOptions: { lang }, globals } });
}
