import * as z from "zod";

import { assetUrl } from "./urls";

const RESULT_LIMIT = 8;

const searchResultDataSchema = z.object({
  url: z.string(),
  excerpt: z.string(),
  meta: z.object({ title: z.string() }).partial(),
  sub_results: z.array(z.object({ title: z.string(), url: z.string(), excerpt: z.string() })).default([]),
});

const searchResponseSchema = z.object({
  results: z.array(
    z.object({ id: z.string(), data: z.function({ input: [], output: z.promise(searchResultDataSchema) }) }),
  ),
});

const pagefindModuleSchema = z.object({
  search: z.function({ input: [z.string()], output: z.promise(searchResponseSchema) }),
});

export type SearchResult = z.infer<typeof searchResultDataSchema>;

export type SearchPages = (query: string) => Promise<SearchResult[]>;

export async function loadPagefind(): Promise<SearchPages> {
  const pagefindUrl = assetUrl("pagefind/pagefind.js");
  const { search } = pagefindModuleSchema.parse(await import(pagefindUrl));

  return async (query) => {
    const { results } = await search(query);

    return Promise.all(results.slice(0, RESULT_LIMIT).map(async (result) => result.data()));
  };
}

export function excerptParts(excerpt: string): ExcerptPart[] {
  return excerpt
    .split(/(<mark>.*?<\/mark>)/u)
    .filter((part) => part !== "")
    .map((part) => {
      const marked = /^<mark>(?<text>.*)<\/mark>$/u.exec(part)?.groups?.["text"];

      return marked === undefined
        ? { text: decodeEntities(part), marked: false }
        : { text: decodeEntities(marked), marked: true };
    });
}

export interface ExcerptPart {
  text: string;
  marked: boolean;
}

const ENTITIES = new Map([
  ["&amp;", "&"],
  ["&lt;", "<"],
  ["&gt;", ">"],
  ["&quot;", '"'],
  ["&#39;", "'"],
  ["&#x27;", "'"],
]);

function decodeEntities(text: string): string {
  return text.replaceAll(/&(?:amp|lt|gt|quot|#39|#x27);/gu, (entity) => ENTITIES.get(entity) ?? entity);
}
