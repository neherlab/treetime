import * as z from "zod";

export const zPickFilesRequest = z.strictObject({
  title: z.string(),
  extensions: z.array(z.string()),
  multiple: z.boolean(),
});

export type PickFilesRequest = z.infer<typeof zPickFilesRequest>;

export const zPickedFiles = z.array(z.string());

export interface LocalFiles {
  pickFiles(request: PickFilesRequest): Promise<string[]>;
  pathForFile(file: File): string;
}
