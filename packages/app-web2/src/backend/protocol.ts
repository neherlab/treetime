import * as z from "zod";

export const BACKEND_WORKER_NAME = "treetime-backend";

export const workerStartSchema = z.strictObject({
  kind: z.literal("start"),
  module: z.instanceof(WebAssembly.Module),
  port: z.instanceof(MessagePort),
  jobs: z.number().int().positive(),
  threadedWasmUrl: z.string().optional(),
});

export const workerEventSchema = z.discriminatedUnion("kind", [
  z.strictObject({ kind: z.literal("ready") }),
  z.strictObject({ kind: z.literal("failed"), message: z.string() }),
  z.strictObject({ kind: z.literal("crashed"), message: z.string() }),
]);

export type WorkerStart = z.infer<typeof workerStartSchema>;

export type WorkerEvent = z.infer<typeof workerEventSchema>;
