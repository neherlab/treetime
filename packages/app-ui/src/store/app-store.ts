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

<<<<<<< HEAD
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
=======
  files: Partial<Record<FileSlotKind, LoadedFile>>;
  setFile: (kind: FileSlotKind, file: LoadedFile | undefined) => void;

  runStatus: RunStatus;
  progress: ProgressInfo | undefined;
  setRunStatus: (status: RunStatus) => void;
  setProgress: (progress: ProgressInfo | undefined) => void;
>>>>>>> 0d74b8c0 (feat(app-ui): add UI mockup with layout shell, input panel, and results placeholders)

  showResults: boolean;
  setShowResults: (show: boolean) => void;

  darkMode: boolean;
  toggleDarkMode: () => void;

  sidebarCollapsed: boolean;
  toggleSidebar: () => void;
}

export const useAppStore = create<AppState>((set) => ({
<<<<<<< HEAD
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
=======
  activeCommand: "timetree",
  setActiveCommand: (command) =>
    set({ activeCommand: command, showResults: false, runStatus: "idle", progress: undefined }),
>>>>>>> 0d74b8c0 (feat(app-ui): add UI mockup with layout shell, input panel, and results placeholders)

  files: {},
  setFile: (kind, file) =>
    set((state) => ({
      files: { ...state.files, [kind]: file },
    })),
<<<<<<< HEAD
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
=======

  runStatus: "idle",
  progress: undefined,
  setRunStatus: (status) => set({ runStatus: status }),
  setProgress: (progress) => set({ progress }),
>>>>>>> 0d74b8c0 (feat(app-ui): add UI mockup with layout shell, input panel, and results placeholders)

  showResults: false,
  setShowResults: (show) => set({ showResults: show }),

  darkMode: false,
  toggleDarkMode: () =>
    set((state) => {
      const next = !state.darkMode;
      document.documentElement.classList.toggle("dark", next);
      return { darkMode: next };
    }),

  sidebarCollapsed: false,
  toggleSidebar: () => set((state) => ({ sidebarCollapsed: !state.sidebarCollapsed })),
}));
