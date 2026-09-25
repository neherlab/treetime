import type { AuspiceThunk } from "auspice/src/state";

export declare function updateVisibleTipsAndBranchThicknesses(options: {
  root: [number | undefined, number | undefined];
}): AuspiceThunk;

export declare function applyFilter(
  mode: "add" | "remove" | "inactivate" | "set",
  trait: string | symbol,
  values: string[],
): AuspiceThunk;
