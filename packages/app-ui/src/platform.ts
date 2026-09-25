import type { LocalFiles } from "@neherlab/app-contracts";
import { createContext, useContext } from "react";

export const LocalFilesContext = createContext<LocalFiles | null>(null);

export function useLocalFiles(): LocalFiles | null {
  return useContext(LocalFilesContext);
}
