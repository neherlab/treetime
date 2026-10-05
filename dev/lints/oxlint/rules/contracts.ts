import type { ESTree } from "@oxlint/plugins";

export const CONTRACTS_PACKAGE = "@neherlab/app-contracts";

export function isContractsSource(source: string): boolean {
  return source === CONTRACTS_PACKAGE || source.startsWith(`${CONTRACTS_PACKAGE}/`);
}

export function contractImports(program: ESTree.Program): Set<string> {
  const names = new Set<string>();

  for (const statement of program.body) {
    if (statement.type !== "ImportDeclaration" || !isContractsSource(statement.source.value)) {
      continue;
    }

    for (const specifier of statement.specifiers) {
      names.add(specifier.local.name);
    }
  }

  return names;
}
