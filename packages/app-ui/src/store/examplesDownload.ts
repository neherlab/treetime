import type { ExamplesDownloadStatus } from "@neherlab/app-contracts";
import { create } from "zustand";

interface ExamplesDownloadState {
  status: ExamplesDownloadStatus | undefined;
  apply: (status: ExamplesDownloadStatus) => void;
}

export const useExamplesDownloadStore = create<ExamplesDownloadState>()((set) => ({
  status: undefined,
  apply: (status) => set((state) => ({ status: newerDownloadStatus(state.status, status) })),
}));

export function newerDownloadStatus(
  current: ExamplesDownloadStatus | undefined,
  incoming: ExamplesDownloadStatus | undefined,
): ExamplesDownloadStatus | undefined {
  if (incoming === undefined) {
    return current;
  }

  if (current === undefined) {
    return incoming;
  }

  return (incoming.seq ?? -1) > (current.seq ?? -1) ? incoming : current;
}
