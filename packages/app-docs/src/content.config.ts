import { glob } from "astro/loaders";
import { defineCollection } from "astro:content";
import * as z from "zod";

const docs = defineCollection({
  loader: glob({ pattern: "**/*.{md,mdx}", base: "./src/content/docs" }),
  schema: z.object({
    title: z.string(),
    description: z.string().optional(),
    template: z.enum(["doc", "splash"]).default("doc"),
    tableOfContents: z.boolean().default(true),
    editUrl: z.literal(false).optional(),
  }),
});

const meta = defineCollection({
  loader: glob({ pattern: "**/meta.json", base: "./src/content/docs" }),
  schema: z.object({
    title: z.string().optional(),
    pages: z.array(z.string()).optional(),
  }),
});

export const collections = { docs, meta };
