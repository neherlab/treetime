import type { JsonValue, ListItemKind } from "@neherlab/app-contracts";

import { isJsonObject } from "./json";
import { parseNumber } from "./numbers";

const TAB_ESCAPE = "\\t";

export function formatList(value: JsonValue | undefined): string {
  return Array.isArray(value) ? value.map((item) => itemText(item).replaceAll("\t", TAB_ESCAPE)).join(" ") : "";
}

export function parseList(text: string, itemKind: ListItemKind): JsonValue[] {
  return text
    .trim()
    .split(/\s+/u)
    .flatMap((token) => {
      if (token === "") {
        return [];
      }

      const item = token.replaceAll(TAB_ESCAPE, "\t");

      return [itemKind === "string" ? item : parseNumber(item)];
    });
}

function itemText(item: JsonValue): string {
  return Array.isArray(item) || isJsonObject(item) ? JSON.stringify(item) : `${item}`;
}
