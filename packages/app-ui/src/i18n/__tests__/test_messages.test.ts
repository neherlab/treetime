import * as z from "zod";

import de from "../../../messages/de.json";
import en from "../../../messages/en.json";
import { m } from "../locale";

const PARAMETER = /\{(?<name>[a-z_]+)\}/gu;

const messageTextsSchema = z.union([
  z.string().transform((text) => [text]),
  z
    .array(z.object({ match: z.record(z.string(), z.string()) }))
    .transform((variants) => variants.flatMap((variant) => Object.values(variant.match))),
]);

const messageFileSchema = z.record(z.string(), messageTextsSchema);

describe("message files", () => {
  test("the German messages have every English key with the same parameters", () => {
    expect(messageParameters(de)).toStrictEqual(messageParameters(en));
  });
});

describe("plural messages", () => {
  test.each([
    ["en", 1, "1 setting matches"],
    ["en", 3, "3 settings match"],
    ["de", 1, "1 Einstellung passt"],
    ["de", 3, "3 Einstellungen passen"],
  ] as const)("the settings filter count in %s for %d", (locale, count, expected) => {
    const text: string = m.settings_filter_count({ count }, { locale });

    expect(text).toBe(expected);
  });
});

function messageParameters(file: unknown): Record<string, string[]> {
  const { $schema: _, ...messages } = messageFileSchema.parse(file);

  return Object.fromEntries(Object.entries(messages).map(([key, texts]) => [key, parameterNames(texts)]));
}

function parameterNames(texts: string[]): string[] {
  return [
    ...new Set(texts.flatMap((text) => [...text.matchAll(PARAMETER)].map((match) => match.groups?.["name"] ?? ""))),
  ].toSorted();
}
