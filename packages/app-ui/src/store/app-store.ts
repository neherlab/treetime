import { create } from "zustand";
import type { CommandName, FileSlotKind, LoadedFile } from "../types";

type RunStatus = "idle" | "running" | "completed" | "failed";

interface ProgressInfo {
  stage: string;
  fraction: number;
  message: string;
}

interface AppState {
  activeCommand: CommandName;
  setActiveCommand: (command: CommandName) => void;

  selectedDataset: string;
  setSelectedDataset: (dataset: string) => void;

  files: Partial<Record<FileSlotKind, LoadedFile>>;
  setFile: (kind: FileSlotKind, file: LoadedFile | undefined) => void;
  resetForm: () => void;

  runStatus: RunStatus;
  progress: ProgressInfo | undefined;
  abortController: AbortController | undefined;
  setRunStatus: (status: RunStatus) => void;
  setProgress: (progress: ProgressInfo | undefined) => void;
  setAbortController: (controller: AbortController | undefined) => void;
  cancelRun: () => void;

  showResults: boolean;
  setShowResults: (show: boolean) => void;

  sidebarCollapsed: boolean;
  toggleSidebar: () => void;
}

export const useAppStore = create<AppState>((set) => ({
  activeCommand: "ancestral",
  setActiveCommand: (command) =>
    set((state) => {
      state.abortController?.abort();
      return {
        activeCommand: command,
        showResults: false,
        runStatus: "idle",
        progress: undefined,
        abortController: undefined,
      };
    }),

  selectedDataset: "",
  setSelectedDataset: (dataset) => set({ selectedDataset: dataset }),

  files: {},
  setFile: (kind, file) =>
    set((state) => ({
      files: { ...state.files, [kind]: file },
    })),
  resetForm: () =>
    set((state) => {
      state.abortController?.abort();
      return {
        files: {},
        selectedDataset: "",
        runStatus: "idle",
        progress: undefined,
        abortController: undefined,
        showResults: false,
      };
    }),

  runStatus: "idle",
  progress: undefined,
  abortController: undefined,
  setRunStatus: (status) => set({ runStatus: status }),
  setProgress: (progress) => set({ progress }),
  setAbortController: (controller) => set({ abortController: controller }),
  cancelRun: () =>
    set((state) => {
      state.abortController?.abort();
      return { runStatus: "idle", progress: undefined, abortController: undefined };
    }),

  showResults: false,
  setShowResults: (show) => set({ showResults: show }),

  sidebarCollapsed: false,
  toggleSidebar: () => set((state) => ({ sidebarCollapsed: !state.sidebarCollapsed })),
}));
