import assert from "node:assert/strict";
import { test } from "node:test";

// oxlint-disable-next-line no-restricted-imports -- the test reads the lint configuration at the repository root, which belongs to no package
import config, { packageGraphPatterns } from "../../../../oxlint.config.ts";

const PROTECTED_RULES = [
  "anti-slop/no-unknown-parameters",
  "anti-slop/no-unknown-returns",
  "anti-slop/no-unsafe-dictionary-type",
  "anti-slop/no-runtime-typeof",
];

const PACKAGE_FILES = /^packages\/(?<dir>[^/]+)\//u;

const overrides = config.overrides ?? [];

void test("no override outside tests and lint code turns off a type-safety rule", () => {
  const offending = overrides
    .filter((override) => !override.files.every(isTestOrLintGlob))
    .flatMap((override) =>
      PROTECTED_RULES.flatMap((rule) =>
        override.rules?.[rule] === "off" ? [`${override.files.join(", ")}: ${rule}`] : [],
      ),
    );

  assert.deepStrictEqual(offending, []);
});

void test("every override that restricts imports keeps the package graph of its package", () => {
  const missing = overrides.flatMap((override) => {
    const rule = override.rules?.["no-restricted-imports"];

    if (rule === undefined) {
      return [];
    }

    const groups = Array.isArray(rule) ? JSON.stringify(rule.slice(1)) : "";

    return override.files.flatMap((files) => {
      const dir = PACKAGE_FILES.exec(files)?.groups?.["dir"];
      const graph = dir === undefined ? [] : packageGraphPatterns(dir);

      return graph.every((pattern) => groups.includes(JSON.stringify(pattern))) ? [] : [files];
    });
  });

  assert.deepStrictEqual(missing, []);
});

void test("the package graph of a UI package bans the packages it does not declare", () => {
  const [graph] = packageGraphPatterns("app-ui");

  assert.deepStrictEqual(
    [graph?.group.includes("@neherlab/app-desktop"), graph?.group.includes("@neherlab/app-contracts")],
    [true, false],
  );
});

function isTestOrLintGlob(glob: string): boolean {
  return glob.startsWith("dev/lints/") || glob.includes("__tests__") || /\.(test|spec)\./u.test(glob);
}
