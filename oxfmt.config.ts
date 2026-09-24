import { defineConfig } from "oxfmt";

export default defineConfig({
  ignorePatterns: [
    "**/*.md",
    "**/*.toml",
    "**/__fixtures__",
    "bun.lock",
    "data",
    "dev/lints/oxlint-anti-slop",
    "kb",
    "packages/app-contracts/openapi.yaml",
    "packages/app-contracts/src/generated",
    "packages/app-napi/index.d.ts",
    "packages/app-output/src/__tests__/schemas",
    "packages/legacy",
    "packages/schemas",
    "test_scripts",
  ],
  printWidth: 120,
  tabWidth: 2,
  useTabs: false,
  semi: true,
  singleQuote: false,
  quoteProps: "as-needed",
  trailingComma: "all",
  sortImports: true,
  sortPackageJson: true,
  sortTailwindcss: {
    functions: ["cn", "clsx", "cva", "tw"],
  },
});
