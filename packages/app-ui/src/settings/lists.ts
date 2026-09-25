import type { ListItemKind } from "@neherlab/app-contracts";

import { isJsonObject, type JsonValue } from "./json";

const TAB_ESCAPE = "\\t";

export function formatList(value: JsonValue): string {
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

      return [itemKind === "string" ? item : Number(item)];
    });
}

function itemText(item: JsonValue): string {
  return Array.isArray(item) || isJsonObject(item) ? JSON.stringify(item) : `${item}`;
}
