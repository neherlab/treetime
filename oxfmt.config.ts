import { defineConfig } from "oxfmt";

export default defineConfig({
  ignorePatterns: [
    "dist",
    "dist-electron",
    "node_modules",
    "bun.lock",
    "*.tsbuildinfo",
    ".turbo",
    ".build",
    "coverage",
    "packages/app-contracts/src/generated",
    "dev/lints/oxlint-anti-slop",
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
