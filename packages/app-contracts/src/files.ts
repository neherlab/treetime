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

export const zPickFolderRequest = z.strictObject({
  title: z.string(),
});

export type PickFolderRequest = z.infer<typeof zPickFolderRequest>;

export const zPickedFolder = z.string().nullable();

export interface WorkspaceShell {
  pickFolder(request: PickFolderRequest): Promise<string | null>;
  restartBackend(): Promise<void>;
}
