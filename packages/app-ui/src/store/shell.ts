import type { AppCommand } from "@neherlab/app-contracts";
import { create } from "zustand";

const COMPARE_LIMIT = 2;

interface ShellState {
  runFilter: string;
  commandFilter: AppCommand | null;
  compareIds: string[];
  paletteOpen: boolean;
  setRunFilter: (runFilter: string) => void;
  setCommandFilter: (commandFilter: AppCommand | null) => void;
  toggleCompare: (id: string) => void;
  clearCompare: () => void;
  setPaletteOpen: (paletteOpen: boolean) => void;
}

export const useShellStore = create<ShellState>()((set) => ({
  runFilter: "",
  commandFilter: null,
  compareIds: [],
  paletteOpen: false,
  setRunFilter: (runFilter) => set({ runFilter }),
  setCommandFilter: (commandFilter) => set({ commandFilter }),
  toggleCompare: (id) =>
    set((state) => ({
      compareIds: state.compareIds.includes(id)
        ? state.compareIds.filter((selected) => selected !== id)
        : [...state.compareIds, id].slice(-COMPARE_LIMIT),
    })),
  clearCompare: () => set({ compareIds: [] }),
  setPaletteOpen: (paletteOpen) => set({ paletteOpen }),
}));
