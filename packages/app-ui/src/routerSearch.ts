import { retainSearchParams } from "@tanstack/react-router";
import * as z from "zod";

const rootSearchSchema = z.object({
  // oxlint-disable-next-line unicorn/prefer-top-level-await -- `catch` of a zod schema replaces an invalid value; it chains no promise
  jobs: z.union([z.number(), z.string()]).optional().catch(undefined),
});

export type RootSearch = z.infer<typeof rootSearchSchema>;

export const rootSearch = {
  validateSearch: rootSearchSchema,
  search: { middlewares: [retainSearchParams<RootSearch>(["jobs"])] },
};
