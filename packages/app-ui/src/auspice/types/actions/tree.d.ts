import type { AuspiceThunk } from "../../state";

export declare function updateVisibleTipsAndBranchThicknesses(options: {
  root: [number | undefined, number | undefined];
}): AuspiceThunk;

export declare function applyFilter(
  mode: "add" | "inactivate" | "remove" | "set" | "focus",
  trait: string | symbol,
  values: readonly string[],
): AuspiceThunk;
