import type { AppCommand } from "@neherlab/app-contracts";
import { isMap, isScalar, parseDocument, stringify } from "yaml";

import { isAppCommand } from "./commands";
import { pathList } from "./inputs";
import { getAt, type JsonObject } from "./json";
import type { SettingSpec } from "./schema";

const SCHEMA_FILE_PATTERN = /input-config-([a-z]+)\.schema\.json/u;

export function configTextCommand(text: string): AppCommand | null {
  const command = SCHEMA_FILE_PATTERN.exec(text)?.[1];

  return isAppCommand(command) ? command : null;
}

export function withMissingInputs(text: string, specs: readonly SettingSpec[], current: JsonObject): string {
  const document = parseDocument(text);

  if (document.errors.length > 0 || !isMap(document.contents)) {
    return text;
  }

  const present = new Set(
    document.contents.items.flatMap((pair) => (isScalar(pair.key) ? [String(pair.key.value)] : [])),
  );

  const missing = specs.flatMap((spec) => {
    const value = getAt(current, spec.path);

    return spec.pathRole === "input" && !present.has(spec.key) && pathList(value).length > 0
      ? [[spec.key, value ?? null] as const]
      : [];
  });

  if (missing.length === 0) {
    return text;
  }

  const separator = text.endsWith("\n") || text === "" ? "" : "\n";

  return `${text}${separator}${stringify(Object.fromEntries(missing))}`;
}
