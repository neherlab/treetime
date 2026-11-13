import * as z from "zod";

const bytesSchema = z.custom<Uint8Array>((value) => value instanceof Uint8Array, "expected bytes");

const blobSchema = z.instanceof(Blob);

export const portHeaderSchema = z.strictObject({ name: z.string(), value: z.string() });

const requestFields = {
  seq: z.number().int().nonnegative(),
  method: z.string(),
  url: z.string(),
  headers: z.array(portHeaderSchema),
};

export const portRequestSchema = z.strictObject({ ...requestFields, body: bytesSchema.exactOptional() });

export const blobPortRequestSchema = z.strictObject({ ...requestFields, body: blobSchema.exactOptional() });

export const portMessageSchema = z.discriminatedUnion("kind", [
  z.strictObject({ kind: z.literal("request"), request: portRequestSchema }),
  z.strictObject({ kind: z.literal("abort"), seq: z.number().int().nonnegative() }),
]);

export const blobPortMessageSchema = z.discriminatedUnion("kind", [
  z.strictObject({ kind: z.literal("request"), request: blobPortRequestSchema }),
  z.strictObject({ kind: z.literal("abort"), seq: z.number().int().nonnegative() }),
]);

export const portReplySchema = z.discriminatedUnion("kind", [
  z.strictObject({
    kind: z.literal("head"),
    seq: z.number().int().nonnegative(),
    status: z.number().int(),
    headers: z.array(portHeaderSchema),
  }),
  z.strictObject({ kind: z.literal("chunk"), seq: z.number().int().nonnegative(), data: bytesSchema }),
  z.strictObject({ kind: z.literal("end"), seq: z.number().int().nonnegative() }),
  z.strictObject({ kind: z.literal("reset"), seq: z.number().int().nonnegative(), message: z.string() }),
]);

export type PortHeader = z.infer<typeof portHeaderSchema>;

export type PortRequest = z.infer<typeof portRequestSchema>;

export type BlobPortRequest = z.infer<typeof blobPortRequestSchema>;

export type AnyPortRequest = PortRequest | BlobPortRequest;

export type PortMessage = z.infer<typeof portMessageSchema>;

export type BlobPortMessage = z.infer<typeof blobPortMessageSchema>;

export type PortReply = z.infer<typeof portReplySchema>;
