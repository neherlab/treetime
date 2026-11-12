import { retainSearchParams } from "@tanstack/react-router";
import * as z from "zod";

const rootSearchSchema = z.object({ jobs: z.union([z.number(), z.string()]).optional() });

export type RootSearch = z.infer<typeof rootSearchSchema>;

export const rootSearch = {
  validateSearch: (search: Record<string, unknown>): RootSearch => {
    const parsed = rootSearchSchema.safeParse(search);

    return parsed.success ? { jobs: parsed.data.jobs } : {};
  },
  search: { middlewares: [retainSearchParams<RootSearch>(["jobs"])] },
};
