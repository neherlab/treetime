import type { SortingState } from "@tanstack/react-table";
import type { SortDescriptor } from "react-aria-components";

export function sortDescriptor(sorting: SortingState): SortDescriptor | undefined {
  const [sort] = sorting;

  if (sort === undefined) {
    return undefined;
  }

  return { column: sort.id, direction: sort.desc ? "descending" : "ascending" };
}
