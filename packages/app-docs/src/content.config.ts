import { glob } from "astro/loaders";
import { defineCollection } from "astro:content";
import * as z from "zod";

import { CONTENT_ROOTS } from "./lib/content-roots";

const pageSchema = z.object({
  title: z.string(),
  description: z.string().optional(),
  template: z.enum(["doc", "splash"]).default("doc"),
  tableOfContents: z.boolean().default(true),
  sourceHash: z.string().optional(),
});

const docs = defineCollection({
  loader: glob({ pattern: "*/**/*.{md,mdx}", base: CONTENT_ROOTS.docs }),
  schema: pageSchema,
});

const generated = defineCollection({
  loader: glob({ pattern: "*/**/*.{md,mdx}", base: CONTENT_ROOTS.generated }),
  schema: pageSchema,
});

const meta = defineCollection({
  loader: glob({ pattern: "*/**/meta.json", base: CONTENT_ROOTS.docs }),
  schema: z.object({
    title: z.string().optional(),
    pages: z.array(z.string()).optional(),
    sourceHash: z.string().optional(),
  }),
});

export const collections = { docs, generated, meta };
