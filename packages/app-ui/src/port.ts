import * as z from "zod";

const bytesSchema = z.custom<Uint8Array>((value) => value instanceof Uint8Array, "expected bytes");

export const portHeaderSchema = z.strictObject({ name: z.string(), value: z.string() });

export const portRequestSchema = z.strictObject({
  seq: z.number().int().nonnegative(),
  method: z.string(),
  url: z.string(),
  headers: z.array(portHeaderSchema),
  body: bytesSchema.exactOptional(),
});

export const portMessageSchema = z.discriminatedUnion("kind", [
  z.strictObject({ kind: z.literal("request"), request: portRequestSchema }),
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

export type PortMessage = z.infer<typeof portMessageSchema>;

export type PortReply = z.infer<typeof portReplySchema>;
