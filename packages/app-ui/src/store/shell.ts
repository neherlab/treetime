import type { AppCommand, JobId } from "@neherlab/app-contracts";
import { create } from "zustand";

interface ShellState {
  runFilter: string;
  commandFilter: AppCommand | null;
  selectedIds: readonly JobId[];
  renamingId: JobId | undefined;
  confirmingDelete: readonly JobId[] | undefined;
  paletteOpen: boolean;
  workspaceOpen: boolean;
  setRunFilter: (runFilter: string) => void;
  setCommandFilter: (commandFilter: AppCommand | null) => void;
  toggleSelected: (id: JobId) => void;
  clearSelection: () => void;
  setRenamingId: (renamingId: JobId | undefined) => void;
  setConfirmingDelete: (confirmingDelete: readonly JobId[] | undefined) => void;
  setPaletteOpen: (paletteOpen: boolean) => void;
  setWorkspaceOpen: (workspaceOpen: boolean) => void;
}

export const useShellStore = create<ShellState>()((set) => ({
  runFilter: "",
  commandFilter: null,
  selectedIds: [],
  renamingId: undefined,
  confirmingDelete: undefined,
  paletteOpen: false,
  workspaceOpen: false,
  setRunFilter: (runFilter) => set({ runFilter }),
  setCommandFilter: (commandFilter) => set({ commandFilter }),
  toggleSelected: (id) =>
    set(({ selectedIds }) => ({
      selectedIds: selectedIds.includes(id) ? selectedIds.filter((selected) => selected !== id) : [...selectedIds, id],
    })),
  clearSelection: () => set({ selectedIds: [] }),
  setRenamingId: (renamingId) => set({ renamingId }),
  setConfirmingDelete: (confirmingDelete) => set({ confirmingDelete }),
  setPaletteOpen: (paletteOpen) => set({ paletteOpen }),
  setWorkspaceOpen: (workspaceOpen) => set({ workspaceOpen }),
}));
