import type { SaveRunArchiveRequest, SaveRunFileRequest } from "@neherlab/app-napi";

export type SaveRunFileDialog = Omit<SaveRunFileRequest, "destination"> & { name: string };

export type SaveRunArchiveDialog = Omit<SaveRunArchiveRequest, "destination"> & { name: string };

export type SaveReply = { saved: boolean } | { error: string };
