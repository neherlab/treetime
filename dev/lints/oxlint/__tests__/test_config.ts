import assert from "node:assert/strict";
import { existsSync, readdirSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { test } from "node:test";

import * as z from "zod";

import config from "../../../../oxlint.config.ts";
import { packageGraphPatterns } from "../config.ts";

const PROTECTED_RULES = [
  "anti-slop/no-unknown-parameters",
  "anti-slop/no-unknown-returns",
  "anti-slop/no-unsafe-dictionary-type",
  "anti-slop/no-runtime-typeof",
];

const PACKAGE_FILES = /^packages\/(?<dir>[^/]+)\//u;

const ROOT = join(import.meta.dirname, "../../../..");

const zManifest = z.object({
  dependencies: z.record(z.string(), z.string()).optional(),
  devDependencies: z.record(z.string(), z.string()).optional(),
  peerDependencies: z.record(z.string(), z.string()).optional(),
});

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
      const graph = dir === undefined ? [] : packageGraphPatterns(ROOT, dir);

      return graph.every((pattern) => groups.includes(JSON.stringify(pattern))) ? [] : [files];
    });
  });

  assert.deepStrictEqual(missing, []);
});

void test("the package graph of a package bans no workspace package it declares", () => {
  const banned = readdirSync(join(ROOT, "packages")).flatMap((dir) => {
    const group = new Set(packageGraphPatterns(ROOT, dir).flatMap((pattern) => pattern.group));

    return declaredWorkspacePackages(dir).filter((name) => group.has(name));
  });

  assert.deepStrictEqual(banned, []);
});

function declaredWorkspacePackages(dir: string): string[] {
  const manifest = join(ROOT, "packages", dir, "package.json");

  if (!existsSync(manifest)) {
    return [];
  }

  const parsed = zManifest.parse(JSON.parse(readFileSync(manifest, "utf8")));

  return Object.keys({ ...parsed.dependencies, ...parsed.devDependencies, ...parsed.peerDependencies }).filter((name) =>
    name.startsWith("@neherlab/"),
  );
}

function isTestOrLintGlob(glob: string): boolean {
  return glob.startsWith("dev/lints/") || glob.includes("__tests__") || /\.(test|spec)\./u.test(glob);
}
