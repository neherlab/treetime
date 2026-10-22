import * as z from "zod";

export const JOB_WORKER_NAME = "treetime-run";

export const CHUNK_WINDOW = 4;

export const moduleVariantSchema = z.enum(["single", "threaded"]);

const bytesSchema = z.custom<Uint8Array>((value) => value instanceof Uint8Array, "expected bytes");

const attempt = z.number().int().positive();

export const hostMessageSchema = z.discriminatedUnion("kind", [
  z.strictObject({
    kind: z.literal("start"),
    attempt,
    variant: moduleVariantSchema,
    module: z.instanceof(WebAssembly.Module),
    job: z.string(),
    jobs: z.number().int().positive(),
  }),
  z.strictObject({ kind: z.literal("input-begin"), attempt, path: z.string() }),
  z.strictObject({ kind: z.literal("input-chunk"), attempt, path: z.string(), data: bytesSchema }),
  z.strictObject({ kind: z.literal("run"), attempt }),
  z.strictObject({ kind: z.literal("output-ack"), attempt }),
]);

export const jobMessageSchema = z.discriminatedUnion("kind", [
  z.strictObject({ kind: z.literal("started"), attempt }),
  z.strictObject({ kind: z.literal("failed"), attempt, lines: z.array(z.string()) }),
  z.strictObject({ kind: z.literal("input-ack"), attempt }),
  z.strictObject({ kind: z.literal("event"), attempt, event: z.string() }),
  z.strictObject({ kind: z.literal("output-begin"), attempt, path: z.string() }),
  z.strictObject({ kind: z.literal("output-chunk"), attempt, path: z.string(), data: bytesSchema }),
  z.strictObject({ kind: z.literal("done"), attempt, outcome: z.string() }),
  z.strictObject({ kind: z.literal("crashed"), attempt, lines: z.array(z.string()) }),
]);

export const messageAttemptSchema = z.object({ attempt });

export type ModuleVariant = z.infer<typeof moduleVariantSchema>;

export type HostMessage = z.infer<typeof hostMessageSchema>;

export type JobMessage = z.infer<typeof jobMessageSchema>;

export function transferList(chunk: Uint8Array): Transferable[] {
  return chunk.buffer instanceof ArrayBuffer ? [chunk.buffer] : [];
}

export class ChunkWindow {
  #free = CHUNK_WINDOW;
  #waiting: Array<(open: boolean) => void> = [];
  #closed = false;

  acquire(): Promise<boolean> {
    if (this.#closed) {
      return Promise.resolve(false);
    }

    if (this.#free > 0) {
      this.#free -= 1;

      return Promise.resolve(true);
    }

    const { promise, resolve } = Promise.withResolvers<boolean>();

    this.#waiting.push(resolve);

    return promise;
  }

  release(): void {
    const next = this.#waiting.shift();

    if (next === undefined) {
      this.#free = Math.min(CHUNK_WINDOW, this.#free + 1);
    } else {
      next(true);
    }
  }

  close(): void {
    this.#closed = true;
    this.#waiting.splice(0).forEach((resolve) => {
      resolve(false);
    });
  }
}
