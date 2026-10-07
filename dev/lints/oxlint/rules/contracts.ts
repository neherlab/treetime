import type { ESTree } from "@oxlint/plugins";
import * as z from "zod";

const zContractsOptions = z.object({ package: z.string() });

export function contractsPackage(options: readonly unknown[]): string | undefined {
  const parsed = zContractsOptions.safeParse(options[0]);

  return parsed.success ? parsed.data.package : undefined;
}

export function isContractsSource(source: string, contracts: string): boolean {
  return source === contracts || source.startsWith(`${contracts}/`);
}

export function contractImports(program: ESTree.Program, contracts: string): Set<string> {
  const names = new Set<string>();

  for (const statement of program.body) {
    if (statement.type !== "ImportDeclaration" || !isContractsSource(statement.source.value, contracts)) {
      continue;
    }

    for (const specifier of statement.specifiers) {
      names.add(specifier.local.name);
    }
  }

  return names;
}
